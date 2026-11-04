"""Deadline-limited, resumable supervisor for the overnight GHZ30 search.

Run this from the code checkout, using a *separate, exclusive* data worktree:

  python search/overnight30_supervisor.py \
    --until 2026-11-05T05:00:00Z --data-repo /workspace/ghz30-overnight \
    --code-root /workspace/threecontext --python /workspace/threecontext-venv/bin/python

One slot explores asymmetric graphs; the others rotate weak symmetries. Only
finished chunk directories and atomic metadata are committed. Active worker
files are never staged. Every new commit is immediately pushed, and a failed
push must succeed before another commit is made. No force pushes are used.

The deadline watchdog is independent of network publication. --no-publish is
intended for local tests; it performs no Git mutations or network operations.
"""
import argparse
from datetime import datetime, timezone
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import socket
import subprocess
import sys
import threading
import time


FAMILIES = ("c2_pairs", "c2_fixed", "c2_mixed", "c3_fixed", "c3_mixed",
            "c5_double", "c5_fixed", "c5_mixed", "c2_onepair", "c2_twopairs", "c3_onecycle")
BRANCH = "experiments/ghz30-20261104"


def utc_now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def parse_time(value):
    result = datetime.fromisoformat(value.replace("Z", "+00:00"))
    if result.tzinfo is None:
        raise ValueError("--until must include a timezone")
    return result.astimezone(timezone.utc).timestamp()


def atomic_json(path, data):
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def load_json(path, default):
    try:
        return json.loads(path.read_text())
    except FileNotFoundError:
        return default


def signal_group(pid, sig):
    try:
        os.killpg(pid, sig)
    except ProcessLookupError:
        pass


class Supervisor:
    def __init__(self, args):
        self.args = args
        self.repo = args.data_repo.resolve()
        self.code = args.code_root.resolve()
        self.root = self.repo / "search" / "overnight30_runs" / args.run_date
        self.root.mkdir(parents=True, exist_ok=True)
        if not self.root.resolve().is_relative_to(self.repo):
            raise ValueError("run directory escapes the data repository")
        self.relative = self.root.relative_to(self.repo)
        self.deadline = parse_time(args.until)
        self.host = socket.gethostname()
        self.lock = open(self.root / ".supervisor.lock", "a+")
        try:
            fcntl.flock(self.lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError("another supervisor holds this run's lock") from None
        (self.root / ".gitignore").write_text(".supervisor.lock\n.control.json\n*.tmp\n")
        self.manifest = load_json(self.root / "manifest.json", {"jobs": [], "created_utc": utc_now()})
        self.control = load_json(self.root / ".control.json", {
            "pending_directories": [], "push_pending": not args.no_publish})
        self.status = load_json(self.root / "status.json", {})
        self.check_other_host()
        self.active = {}
        self.active_lock = threading.Lock()
        self.stopping = threading.Event()
        self.finished = threading.Event()
        self.reason = None
        self.next_id = max((j["id"] for j in self.manifest["jobs"]), default=-1) + 1
        self.last_publish_attempt = 0.0
        self.last_heartbeat = 0.0
        self.status.update({"host": self.host, "supervisor_pid": os.getpid(),
                            "started_utc": utc_now(), "deadline_utc": args.until,
                            "stopped": False, "stop_reason": None,
                            "branch": args.branch, "publishing_enabled": not args.no_publish})
        if not args.no_publish:
            self.check_repository()
        self.recover_interrupted()
        self.save()

    def check_other_host(self):
        """A fresh clone must not duplicate a recently checkpointed remote run."""
        if self.status.get("stopped") or self.status.get("host") in (None, self.host):
            return
        heartbeat = self.status.get("heartbeat_utc")
        if heartbeat and time.time()-parse_time(heartbeat) < self.args.stale_host_seconds:
            raise RuntimeError("a recent checkpoint reports a supervisor on another host; inspect it first")

    def git(self, *arguments, timeout=30):
        return subprocess.run(["git", *arguments], cwd=self.repo, text=True,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              timeout=timeout, check=False)

    def check_repository(self):
        top = self.git("rev-parse", "--show-toplevel")
        if top.returncode or Path(top.stdout.strip()).resolve() != self.repo:
            raise ValueError("--data-repo must be a Git worktree root")
        branch = self.git("symbolic-ref", "--short", "HEAD")
        if branch.returncode or branch.stdout.strip() != self.args.branch:
            raise ValueError("data worktree is not on the expected dedicated branch")
        remote = self.git("remote", "get-url", "--push", "origin")
        if remote.returncode or remote.stdout.strip() not in (
                "https://github.com/ManekiMeow/threecontext.git",
                "https://github.com/ManekiMeow/threecontext",
                "git@github.com:ManekiMeow/threecontext.git"):
            raise ValueError("unexpected data-worktree push destination")
        staged = self.git("diff", "--cached", "--name-only")
        if staged.returncode or any(not p.startswith(str(self.relative)+"/")
                                    for p in staged.stdout.splitlines()):
            raise RuntimeError("unrelated staged paths exist in the data worktree")

    def recover_interrupted(self):
        for job in self.manifest["jobs"]:
            if job.get("state") != "running":
                continue
            pid = job.get("pid")
            if job.get("host") == self.host and pid:
                try:
                    command = Path(f"/proc/{pid}/cmdline").read_bytes().replace(b"\0", b" ").decode()
                except (OSError, UnicodeError):
                    command = ""
                if "overnight30_" in command and str(self.repo / job["directory"]) in command:
                    signal_group(pid, signal.SIGTERM)
                    # Avoid archiving files while an orphan still writes them.
                    for _ in range(20):
                        if not Path(f"/proc/{pid}").exists():
                            break
                        time.sleep(0.1)
                    signal_group(pid, signal.SIGKILL)
            job.update(state="interrupted", finished_utc=utc_now(), returncode=None,
                       termination_reason="supervisor_restart")
            self.archive(job)
            self.queue_directory(job["directory"])

    def save(self):
        with self.active_lock:
            active = [{"id": item["job"]["id"], "pid": item["process"].pid,
                       "family": item["job"].get("family"), "kind": item["job"]["kind"],
                       "directory": item["job"]["directory"]} for item in self.active.values()]
        self.status.update(heartbeat_utc=utc_now(), active=active,
                           jobs_started=len(self.manifest["jobs"]),
                           jobs_finished=sum(j.get("state") != "running" for j in self.manifest["jobs"]))
        atomic_json(self.root / "manifest.json", self.manifest)
        atomic_json(self.root / ".control.json", self.control)
        atomic_json(self.root / "status.json", self.status)
        self.last_heartbeat = time.monotonic()

    def queue_directory(self, directory):
        if directory not in self.control["pending_directories"]:
            self.control["pending_directories"].append(directory)

    def archive(self, job):
        directory = self.repo / job["directory"]
        for path in directory.rglob("*.jsonl"):
            destination = path.with_suffix(path.suffix + ".gz")
            temporary = destination.with_name(destination.name + ".tmp")
            with open(path, "rb") as source, open(temporary, "wb") as raw:
                with gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as compressed:
                    shutil.copyfileobj(source, compressed)
            temporary.replace(destination)
            path.unlink()
        job["archive_bytes"] = sum(p.stat().st_size for p in directory.rglob("*") if p.is_file())

    def push(self):
        pushed = self.git("push", "origin", f"HEAD:refs/heads/{self.args.branch}")
        if pushed.returncode:
            self.status["publication_error"] = pushed.stdout[-1500:]
            return False
        self.control["push_pending"] = False
        self.status["publication_error"] = None
        self.status["last_pushed_commit"] = self.git("rev-parse", "HEAD").stdout.strip()
        self.status["last_published_utc"] = utc_now()
        return True

    def publish(self, final=False):
        if self.args.no_publish:
            self.control["pending_directories"].clear()
            self.save()
            return True
        self.last_publish_attempt = time.monotonic()
        try:
            # A previous commit must be delivered before another is created.
            if self.control["push_pending"] and not self.push():
                self.save()
                return False
            pending = list(self.control["pending_directories"])
            if not pending and not final:
                self.save()
                return True
            self.save()
            paths = pending + [str(self.relative / name) for name in
                               ("manifest.json", "status.json", ".gitignore")]
            # The dedicated worktree's index is reserved for this supervisor.
            self.check_repository()
            added = self.git("add", "--", *paths)
            if added.returncode:
                raise RuntimeError(added.stdout)
            changed = self.git("diff", "--cached", "--quiet")
            if changed.returncode == 1:
                self.control["push_pending"] = True
                self.save()  # persist recovery intent before creating the commit
                message = "Finish" if final else "Checkpoint"
                committed = self.git("commit", "-m", f"{message} overnight GHZ30 search: "
                                     f"{self.status['jobs_finished']} chunks")
                if committed.returncode:
                    self.control["push_pending"] = False
                    raise RuntimeError(committed.stdout)
                self.control["pending_directories"].clear()
                if not self.push():
                    self.save()
                    return False
            elif changed.returncode:
                raise RuntimeError(changed.stdout)
            else:
                self.control["pending_directories"].clear()
            self.save()
            return True
        except (subprocess.TimeoutExpired, RuntimeError) as error:
            self.status["publication_error"] = str(error)[-1500:]
            self.save()
            return False

    def previous_job(self, kind, family=None):
        for job in reversed(self.manifest["jobs"]):
            if job["kind"] == kind and job.get("family") == family and job.get("state") != "running":
                prior = self.repo / job["directory"]
                checkpoint = (prior / family / f"seed_{job['seed']}" / "state.json"
                              if family else prior / f"asymmetric_seed{job['seed']}.json")
                if checkpoint.exists():
                    return job
        return None

    def next_family(self):
        with self.active_lock:
            occupied = {item["job"].get("family") for item in self.active.values()}
        choices = [f for f in self.args.families if f not in occupied]
        if not choices:
            return None
        return min(choices, key=lambda f: sum(j.get("family") == f for j in self.manifest["jobs"]))

    def launch(self, slot):
        remaining = self.deadline-time.time()
        if remaining < 5:
            return False
        kind = "asymmetric" if slot == 0 else "symmetry"
        family = self.next_family() if kind == "symmetry" else None
        if kind == "symmetry" and family is None:
            return False
        seed = self.args.seed + (self.args.families.index(family)+1 if family else 0)
        identifier = self.next_id
        self.next_id += 1
        label = family or kind
        directory = self.root / "chunks" / f"{identifier:05d}_{label}"
        directory.mkdir(parents=True)
        chunk_seconds = (self.args.first_chunk_seconds if identifier < self.args.workers
                         and self.args.first_chunk_seconds else self.args.chunk_seconds)
        budget = min(float(chunk_seconds), max(1.0, remaining-2))
        command = [str(self.args.python), str(self.code / "search" / f"overnight30_{kind}.py"),
                   "--seconds", str(budget), "--seed", str(seed), "--output-dir", str(directory)]
        if family:
            command += ["--family", family]
        previous = self.previous_job(kind, family)
        if previous:
            prior = self.repo / previous["directory"]
            if kind == "symmetry":
                command += ["--resume-from", str(prior)]
            else:
                checkpoint = prior / f"asymmetric_seed{seed}.json"
                if checkpoint.exists():
                    command += ["--resume-from", str(checkpoint)]
        environment = os.environ.copy()
        for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
            environment[name] = "1"
        environment["PYTHONUNBUFFERED"] = "1"
        log = open(directory / "worker.log", "ab", buffering=0)
        process = subprocess.Popen(command, cwd=self.code, env=environment, stdin=subprocess.DEVNULL,
                                   stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        job = {"id": identifier, "kind": kind, "family": family, "seed": seed,
               "directory": str(directory.relative_to(self.repo)), "started_utc": utc_now(),
               "host": self.host, "pid": process.pid, "seconds_budget": budget, "state": "running",
               "resume_job": previous["id"] if previous else None}
        code_head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=self.code,
                                   capture_output=True, text=True, timeout=5, check=False)
        job["code_commit"] = code_head.stdout.strip() if code_head.returncode == 0 else None
        job["source_sha256"] = {str(path.relative_to(self.code)):
                              hashlib.sha256(path.read_bytes()).hexdigest() for path in
                              (self.code / "search" / f"overnight30_{kind}.py",
                               self.code / "search" / "overnight30_supervisor.py") if path.exists()}
        self.manifest["jobs"].append(job)
        with self.active_lock:
            self.active[slot] = {"job": job, "process": process, "log": log,
                                 "hard_stop": min(self.deadline, time.time()+budget+15),
                                 "term_sent": None, "started_monotonic": time.monotonic()}
        self.save()
        return True

    def watchdog(self):
        while not self.finished.wait(0.25):
            now = time.time()
            if now >= self.deadline:
                self.reason = self.reason or "deadline"
                self.stopping.set()
            with self.active_lock:
                for item in self.active.values():
                    if item["process"].poll() is not None:
                        continue
                    if self.stopping.is_set() or now >= item["hard_stop"]:
                        if item["term_sent"] is None:
                            item["term_sent"] = now
                            item["job"]["termination_reason"] = self.reason or "chunk_timeout"
                            signal_group(item["process"].pid, signal.SIGTERM)
                        elif now-item["term_sent"] > 3:
                            signal_group(item["process"].pid, signal.SIGKILL)

    def reap(self):
        completed = []
        with self.active_lock:
            for slot, item in list(self.active.items()):
                code = item["process"].poll()
                if code is not None:
                    completed.append(item)
                    del self.active[slot]
                    item["job"].update(state="completed" if code == 0 else "failed",
                                       finished_utc=utc_now(), returncode=code,
                                       elapsed_seconds=time.monotonic()-item["started_monotonic"])
        for item in completed:
            # The group can outlive its leader (for example a numerical child
            # ignores SIGTERM). Do not leave those descendants after a batch.
            signal_group(item["process"].pid, signal.SIGKILL)
            item["log"].close()
            self.archive(item["job"])
            self.queue_directory(item["job"]["directory"])
            recent = [job for job in reversed(self.manifest["jobs"])
                      if job["kind"] == item["job"]["kind"] and job["state"] != "running"][:3]
            if len(recent) == 3 and all(job["state"] == "failed" and
                    job.get("elapsed_seconds", 100) < 10 and not job.get("termination_reason")
                    for job in recent):
                self.reason = "repeated_worker_startup_failure"
                self.stopping.set()
        if completed:
            self.save()
        return completed

    def stop(self, signum=None, frame=None):
        self.reason = self.reason or (f"signal_{signum}" if signum else "requested")
        self.stopping.set()

    def run(self):
        signal.signal(signal.SIGTERM, self.stop)
        signal.signal(signal.SIGINT, self.stop)
        watcher = threading.Thread(target=self.watchdog, name="deadline-watchdog", daemon=True)
        watcher.start()
        try:
            while True:
                self.reap()
                if time.time() >= self.deadline:
                    self.reason = self.reason or "deadline"
                    self.stopping.set()
                if not self.stopping.is_set():
                    for slot in range(self.args.workers):
                        with self.active_lock:
                            free = slot not in self.active
                        if free:
                            self.launch(slot)
                if self.control["pending_directories"] and (
                        not self.status.get("publication_error") or
                        time.monotonic()-self.last_publish_attempt >= 60):
                    self.publish()
                elif self.control["push_pending"] and time.monotonic()-self.last_publish_attempt >= 60:
                    self.publish()
                if time.monotonic()-self.last_heartbeat >= 10:
                    self.save()
                archive_bytes = sum(j.get("archive_bytes", 0) for j in self.manifest["jobs"])
                if archive_bytes > self.args.max_archive_mb*1024*1024:
                    self.reason = "archive_size_limit"
                    self.stopping.set()
                with self.active_lock:
                    active = bool(self.active)
                if self.stopping.is_set() and not active:
                    break
                time.sleep(0.5)
            self.status.update(stopped=True, stop_reason=self.reason, finished_utc=utc_now())
            self.save()
            self.publish(final=True)
            print(json.dumps(self.status, indent=2), flush=True)
            return 1 if self.status.get("publication_error") else 0
        finally:
            self.stopping.set()
            with self.active_lock:
                for item in self.active.values():
                    signal_group(item["process"].pid, signal.SIGKILL)
                    item["log"].close()
            self.finished.set()
            watcher.join(timeout=2)
            fcntl.flock(self.lock.fileno(), fcntl.LOCK_UN)
            self.lock.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--until", required=True)
    parser.add_argument("--data-repo", type=Path, required=True)
    parser.add_argument("--code-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--python", type=Path, default=Path(sys.executable))
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--chunk-seconds", type=float, default=600)
    parser.add_argument("--first-chunk-seconds", type=float,
                        help="shorter first batch to verify publication promptly")
    parser.add_argument("--run-date", default="2026-11-04")
    parser.add_argument("--seed", type=int, default=20261104)
    parser.add_argument("--families", type=lambda x: tuple(x.split(",")), default=FAMILIES)
    parser.add_argument("--branch", default=BRANCH)
    parser.add_argument("--stale-host-seconds", type=float, default=1200)
    parser.add_argument("--max-archive-mb", type=float, default=50)
    parser.add_argument("--no-publish", action="store_true")
    args = parser.parse_args()
    if args.workers < 1 or args.workers > 4 or args.chunk_seconds <= 0:
        parser.error("workers must be 1..4 and chunk-seconds must be positive")
    if args.first_chunk_seconds is not None and args.first_chunk_seconds <= 0:
        parser.error("first-chunk-seconds must be positive")
    if not args.families or any(not f for f in args.families) or len(set(args.families)) != len(args.families):
        parser.error("families must be distinct nonempty names")
    if not args.run_date or "/" in args.run_date or args.run_date in (".", ".."):
        parser.error("run-date must be a directory basename")
    sys.exit(Supervisor(args).run())


if __name__ == "__main__":
    main()
