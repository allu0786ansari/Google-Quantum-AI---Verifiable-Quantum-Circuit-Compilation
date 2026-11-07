// Bit-sliced value simulator for kickmix kmx ops; per-Toffoli firing statistics.
// usage: sim ops.bin nops nq nb W init_q.bin init_b.bin seed out_stats.bin
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef uint64_t u64;
static u64 rng_s[2];
static u64 next(void){u64 s1=rng_s[0];const u64 s0=rng_s[1];rng_s[0]=s0;s1^=s1<<23;rng_s[1]=s1^s0^(s1>>17)^(s0>>26);return rng_s[1]+s0;}
int main(int argc,char**argv){
  if(argc<10){fprintf(stderr,"args\n");return 1;}
  long nops=atol(argv[2]); int nq=atoi(argv[3]); long nb=atol(argv[4]); int W=atoi(argv[5]);
  int32_t*ops=malloc(sizeof(int32_t)*5*nops); FILE*f=fopen(argv[1],"rb"); fread(ops,4,5*nops,f); fclose(f);
  u64*Q=calloc((size_t)nq*W,8),*B=calloc((size_t)nb*W,8);
  f=fopen(argv[6],"rb"); fread(Q,8,(size_t)nq*W,f); fclose(f);
  f=fopen(argv[7],"rb"); fread(B,8,(size_t)nb*W,f); fclose(f);
  rng_s[0]=strtoull(argv[8],0,10)|1; rng_s[1]=0x9E3779B97F4A7C15ULL;
  // stats per op: active, a1, b1, fire  (only for CCX/CCZ/CSWAP)
  uint32_t*st=calloc((size_t)nops*4,4);
  int stack[256]; int sp=0; u64*act=malloc(8*W); u64*tmp=malloc(8*W);
  for(long i=0;i<nops;i++){
    int32_t*o=ops+5*i; int op=o[0],a=o[1],b=o[2],c=o[3],cd=o[4];
    if(op==11){stack[sp++]=cd;continue;} if(op==12){sp--;continue;}
    if(op==0||op==13||op==14||op==16) { if(op!=14 && op!=13) continue; }
    for(int w=0;w<W;w++){u64 m=~0ULL; for(int s=0;s<sp;s++) m&=B[(size_t)stack[s]*W+w]; if(cd>=0) m&=B[(size_t)cd*W+w]; act[w]=m;}
    u64*qa=a>=0?Q+(size_t)a*W:0,*qb=b>=0?Q+(size_t)b*W:0,*qc=c>=0?Q+(size_t)c*W:0;
    switch(op){
      case 1: for(int w=0;w<W;w++) qb[w]^=act[w]&qa[w]; break;
      case 2: { uint32_t ac=0,a1=0,b1=0,fi=0; for(int w=0;w<W;w++){u64 m=act[w]; ac+=__builtin_popcountll(m); a1+=__builtin_popcountll(m&qa[w]); b1+=__builtin_popcountll(m&qb[w]); u64 x=m&qa[w]&qb[w]; fi+=__builtin_popcountll(x); qc[w]^=x;} uint32_t*s=st+4*i; s[0]+=ac;s[1]+=a1;s[2]+=b1;s[3]+=fi; } break;
      case 15: { uint32_t ac=0,a1=0,b1=0,fi=0; for(int w=0;w<W;w++){u64 m=act[w]; ac+=__builtin_popcountll(m); a1+=__builtin_popcountll(m&qa[w]); b1+=__builtin_popcountll(m&qb[w]); fi+=__builtin_popcountll(m&qa[w]&qb[w]&qc[w]);} uint32_t*s=st+4*i; s[0]+=ac;s[1]+=a1;s[2]+=b1;s[3]+=fi; } break;
      case 3: for(int w=0;w<W;w++) qa[w]^=act[w]; break;
      case 4: for(int w=0;w<W;w++){u64 d=(qa[w]^qb[w])&act[w]; qa[w]^=d; qb[w]^=d;} break;
      case 5: { uint32_t ac=0,c1=0; for(int w=0;w<W;w++){u64 m=act[w]&qa[w]; ac+=__builtin_popcountll(act[w]); c1+=__builtin_popcountll(m); u64 d=(qb[w]^qc[w])&m; qb[w]^=d; qc[w]^=d;} uint32_t*s=st+4*i; s[0]+=ac;s[1]+=c1;s[3]+=c1; } break;
      case 6: for(int w=0;w<W;w++) qa[w]&=~act[w]; break;
      case 7: { u64*bb=B+(size_t)b*W; for(int w=0;w<W;w++){u64 r=next(); bb[w]=(bb[w]&~act[w])|(r&act[w]); qa[w]&=~act[w];} } break;
      case 8: { u64*bb=B+(size_t)a*W; for(int w=0;w<W;w++) bb[w]^=act[w]; } break;
      case 9: { u64*bb=B+(size_t)a*W; for(int w=0;w<W;w++) bb[w]&=~act[w]; } break;
      case 10:{ u64*bb=B+(size_t)a*W; for(int w=0;w<W;w++) bb[w]|=act[w]; } break;
      default: break;
    }
  }
  f=fopen(argv[9],"wb"); fwrite(st,4,(size_t)nops*4,f); fclose(f);
  // final check: nonzero ancillas (q>=840)
  long bad=0; for(int q=840;q<nq;q++) for(int w=0;w<W;w++) bad+=__builtin_popcountll(Q[(size_t)q*W+w]);
  fprintf(stderr,"done; nonzero ancilla bits at end: %ld\n",bad);
  // dump final syndromes for check
  f=fopen("final_q.bin","wb"); fwrite(Q,8,(size_t)840*W,f); fclose(f);
  return 0;}
