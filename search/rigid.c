/* Filter graph6 input: keep graphs that are 3-colourable and "3-rigid":
   in EVERY proper 3-colouring, every vertex has neighbours of both other colours.
   (Necessary for a vertex-minimal 3-context GHZ graph, see notes.) */
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#define MAXN 64
typedef uint64_t set;
int n; set adj[MAXN]; int col[MAXN]; int ncol; long found;
int ok_flag, any_col;
static int parse(char *s){
  if(s[0]=='>' ) s+=10;
  n=s[0]-63; int k=1;
  for(int i=0;i<n;i++) adj[i]=0;
  int bit=0; int val=0; int p=1;
  for(int j=1;j<n;j++) for(int i=0;i<j;i++){
    if(bit==0){ val=s[p++]-63; bit=6; }
    bit--; if((val>>bit)&1){ adj[i]|=1ULL<<j; adj[j]|=1ULL<<i; }
  }
  return n;
}
static int check_rigid(void){
  for(int v=0;v<n;v++){ int seen=0;
    for(int w=0;w<n;w++) if(adj[v]>>w&1) seen|=1<<col[w];
    seen &= ~(1<<col[v]);
    if(seen!=(7&~(1<<col[v]))) return 0; }
  return 1;
}
/* returns 1 to abort (non-rigid colouring found) */
static int bt(int v,int maxc){
  if(v==n){ any_col=1; if(!check_rigid()) return 1; return 0; }
  for(int c=0;c<3 && c<=maxc+1;c++){
    int okc=1; set a=adj[v];
    while(a){ int w=__builtin_ctzll(a); a&=a-1; if(w<v && col[w]==c){okc=0;break;} }
    if(!okc) continue;
    col[v]=c; if(bt(v+1, c>maxc?c:maxc)) return 1;
  }
  return 0;
}
int main(){
  char line[4096];
  while(fgets(line,sizeof line,stdin)){
    line[strcspn(line,"\n")]=0; parse(line);
    any_col=0;
    if(bt(0,-1)) continue;
    if(!any_col) continue;
    puts(line);
  }
}
