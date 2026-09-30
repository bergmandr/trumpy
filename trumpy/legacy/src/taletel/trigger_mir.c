/* From westerhoff@nevis1.nevis.columbia.eduTue Nov 21 13:52:02 2000

What we have now is
a clustering program that sorts all hits into clusters of adjacent
hits. After that is done, the trigger requirement for a mirror is at least
one cluster with more than 3 hits. 
For clusters with between 3 and 10 hits, we also require that the cluster
has a non-zero time spread sigma_t.
The programs are of course assembler code for the dsp, but for my own
testing, I have a C version which I will attach to this mail.
It reads in the event form a list that looks like this

7 5                mirror 5, 7 hits
2 0 1 54           first hit: slave, channel, pes, time
3 3 1 265          second hit:    dto
4 1 1 37
5 7 1 41
6 9 1 38
9 1 1 276
9 2 1 70           seventh hit: dto      */

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <string.h>
/* caozh Arguments  caozh */
struct arg{
  int nhit;
  int flag[256];
  int slave[256];
  int bdr_chn[256];
  int npe[256];
  int tav[256];
}to_mirtrig_;

int trigger_mir_( void )
{ 
//   printf("now in trigger_mir_\n");
  int i,j,k;
  int hits,nhit;
  int nhis;
  int r3,rr3,r31,r4,r5;
  int nhitscan;
  int nhit_tot;
  int nEvents=0,nTriggers=0,nReject=0;
  int nClusters=0,rClusters=0;
  int nClustersPrint=0;
  int nPrint=50;
  int nclus,rclus;
  int success,odd;
  int start;
  int etrigger,ctrigger;
//   int mirror;

  float mx,my,mx2,my2,mxy;
  float sigma_x2,sigma_y2,sigma_xy;
  float nsigma_x2,nsigma_y2,nsigma_xy;
  float mt,mt2;
  float sigma_t2,sigma_t;
  float nsigma_t2;
  float x,y,d,nd;
  float /*slave,channel,*/time;
  float pes,pes_tot;
  float delta_slave,delta_channel;
  float s,sum;
  float ns,nsum;
  float s_norm,ss;
  float length,width,ratio;
  float pes_cut=0.1;

  typedef struct {
     int flag;
     float slave,channel;
     float time,pes;
  } EVENT;

  EVENT hitbuff[256];
  EVENT clusters[256];
  int clist[48];

  nhit_tot = to_mirtrig_.nhit;

  memset(hitbuff,0,sizeof(hitbuff));
  memset(clusters,0,sizeof(clusters));
  memset(clist,0,sizeof(clist));

  nEvents++;
  nclus = rclus = 0;
  nhitscan = 0;
//   printf("nhit_tot %d\n",nhit_tot);
  for (k=0; k<nhit_tot; k++) {
//     printf("to_mirtrig_.npe[%d] = %d, pes_cut %f\n",k,to_mirtrig_.npe[k],pes_cut);
    if (to_mirtrig_.npe[k] > pes_cut) {
      hitbuff[nhitscan].slave   = to_mirtrig_.slave[k];
      hitbuff[nhitscan].channel = to_mirtrig_.bdr_chn[k];
      hitbuff[nhitscan].time    = to_mirtrig_.tav[k];
      hitbuff[nhitscan].pes     = to_mirtrig_.npe[k];   
      hitbuff[nhitscan].flag    = 0;
      nhitscan++;
    }
  }
/*
 *  sort by cluster
 */

  r31 = r4 = r5 = 0;
//  printf("nhitscan %d\n",nhitscan);
  for (r3=0; r3<nhitscan; r3++) {
    if (!hitbuff[r3].flag) {
      r31 = r3;
      hits = 1;
      hitbuff[r3].flag     = 1;
      clusters[r4].flag    = nclus;
      clusters[r4].slave   = hitbuff[r3].slave;
      clusters[r4].channel = hitbuff[r3].channel;
      clusters[r4].time    = hitbuff[r3].time;
      clusters[r4].pes     = hitbuff[r3].pes;
      r4++;
      
      while (r4 > r5) {
	for (rr3=r31; rr3<nhitscan; rr3++) {
	  if (!hitbuff[rr3].flag) {
	    
	    delta_slave   = hitbuff[rr3].slave - clusters[r5].slave;
	    delta_channel = fabs(hitbuff[rr3].channel - clusters[r5].channel);

	    odd = 0;
	    odd = (((int)clusters[r5].channel)%2);
	    
	    success = 0;
	    if (delta_slave == 0. && delta_channel <= 2) success = 1;
	    if (fabs(delta_slave) == 1. && delta_channel <= 1) success = 1;
	    if (delta_channel == 1.) {
	      if (odd == 1 && delta_slave == -2.) success = 1;
	      if (odd == 0 && delta_slave ==  2.) success = 1;
	    }

	    if (success) {
	      hitbuff[rr3].flag    = 1;
	      clusters[r4].flag    = nclus;
	      clusters[r4].slave   = hitbuff[rr3].slave;
	      clusters[r4].channel = hitbuff[rr3].channel;
	      clusters[r4].time    = hitbuff[rr3].time;
	      clusters[r4].pes     = hitbuff[rr3].pes;
	      r4++;
	      hits++;
	    }
	  }
	}
	r5++;
      }     
      clist[nclus] = hits;
      nclus++;
      nClusters++;
    }
  }       
/*
 *  calculate moments
 */

  etrigger = 0;
  start = 0;
//    printf("nclus %d\n",nclus);
  for (i=0; i<nclus; i++) {
    mx = my = mxy = mt = 0.;
    mx2 = my2 = mt2 = 0.;
    pes_tot = 0.;
    nhit = 0;
    ctrigger = 0;
    
    for (j=start; j<(start+clist[i]); j++) {
      pes  = clusters[j].pes;
      time = clusters[j].time;
      x = 1.+clusters[j].slave;
      y = 1.+clusters[j].channel;
      if (pes > 0.00001) {
	mx += pes * x;
	my += pes * y;
	mt += pes * time;
	mx2 += pes * x * x;
	my2 += pes * y * y;
	mxy += pes * x * y;
	mt2 += pes * time * time;
	pes_tot += pes;
	nhit++;
      }
    }
    
    /* if (clist[i] >= 3) {  this was wrong according to John B. */
     if (clist[i] > 3) {
      rClusters++;
      rclus++;
      
      nsigma_x2 = pes_tot * mx2 - mx * mx;
      nsigma_y2 = pes_tot * my2 - my * my;
      /* debugging AZ */
//       printf("pes_tot %f,mx2 %f,my2 %f,mx %f,my %f\n",
// 	     pes_tot,mx2,my2,mx,my);    
//       printf("nsigma_x2 %f,nsigma_y2 %f\n",
// 	     nsigma_x2,nsigma_y2);

      nsigma_xy = pes_tot * mxy - mx * my;
      nsigma_t2 = pes_tot * mt2 - mt * mt;
      
      sigma_x2 = nsigma_x2/( pes_tot * pes_tot ); 
      /* sigma_x2 = pow(0.9,2.);  pow(0.9,2) doesn't work in Rutgers! */  
      sigma_y2 = nsigma_y2/( pes_tot * pes_tot );
      sigma_xy = nsigma_xy/( pes_tot * pes_tot );
      sigma_t2 = nsigma_t2/( pes_tot * pes_tot );
      sigma_t  = sqrt(sigma_t2);
/*      printf("sigma_x2 %f,sigma_y2 %f\n",
	     sigma_x2,sigma_y2)*/;

      nd   = nsigma_y2 - nsigma_x2;
      nsum = nsigma_y2 + nsigma_x2;
      d    = sigma_y2 - sigma_x2;
      sum  = sigma_y2 + sigma_x2;
      
      ns = nd * nd + 4. * nsigma_xy * nsigma_xy;
      s  = d * d + 4. * sigma_xy * sigma_xy;
      ss  = sqrt(s);
      /* try to prevent crash, AZ */
      if (sum == 0.00000) {
         sum = sum + 0.00001;
         printf("$$$$$$$ sum=0 in trigger_mir $$$$\n");
	   }
      s_norm = ss/sum;
      length = sqrt(0.5*(sum + ss));
      width  = sqrt(0.5*(sum - ss));
      if (width > 0.001) { 
	ratio = length/width;
      } else {
	ratio = 1000.;
      }
      
      if (clist[i] < 10) {
	if (sigma_t > 5.) etrigger = 1; 
      } 
      else {
	etrigger = 1;
	/*            if (s_norm > 0.4) etrigger = 1;    */
      }
      /*
      printf("--------------\n");
      printf("mirror  %d:\n",mirror);
      printf("cluster %d: %d hits\n",i,clist[i]);
      printf("total pes = %f\n",pes_tot);
      printf("s_x + s_y = %f (%f)\n",sum,pow(nsum,2));
      printf("sigma_t   = %f\n",sigma_t);
      printf("sigma_t2  = %f (%f)\n",sigma_t2,nsigma_t2);
      printf("s^2       = %f (%f)\n",s,ns);
      printf("s_norm    = %f\n",s_norm);
      printf("ratio     = %f (%f/%f)\n",ratio,length,width);
      */
      if (clist[i] >= 10 && sigma_t > 5.) ctrigger=1;
      
      if (ctrigger) {
	if (nClustersPrint < nPrint) {
	  nClustersPrint++;
	  nhis = 200-1+nClustersPrint;
	  /*printf("%d %f %f\n",nhis,ratio,s);*/
	}
      }
    }
    start += clist[i];
  }

  /*for (i=0; i<nhitscan; i++) {
    printf("%d %d / (%d) %d %d %f %f \n",(int)hitbuff[i].slave,
	   (int)hitbuff[i].channel,
	   clusters[i].flag,
	   (int)clusters[i].slave,
	   (int)clusters[i].channel,
	   clusters[i].pes,clusters[i].time);
	   }*/
  if (etrigger) {
    nTriggers++;
  } 
  else { 
    nReject++;
    /*if (nReject < nPrint) {
      printf("-----------------rejected--- \n");
      }*/
  }
  /*
  printf("total number of events  : %d\n",nEvents);
  printf("             of triggers: %d (%d)\n",nTriggers,nReject);
  printf("             of clusters: %d (%d)\n",rClusters,nClusters);
  */
//   printf("end of trigger_mir_: nTriggers = %d\n",nTriggers);
  return nTriggers; 
}









