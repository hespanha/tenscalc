/*
  This file is part of Tencalc.

  Copyright (C) 2012-21 The Regents of the University of California
  (author: Dr. Joao Hespanha).  All rights reserved.
*/

void profilingViewFP(FILE *fp)
{
  int nGroups=sizeof(countCallGroup)/sizeof(countCallGroup[0]);
  int nGets  =sizeof(getNames)/sizeof(getNames[0]);
  int nSets  =sizeof(setNames)/sizeof(setNames[0]);
  int nCopies=sizeof(copyNames)/sizeof(copyNames[0]);
  int nFlops =sizeof(flopsNames)/sizeof(flopsNames[0]);

  fprintf(fp,"Operation    Count\n");
  int64_t c=0;
  for (int i=0;i<nFlops;i++) {
    fprintf(fp,"%-10s%9"PRId64"\n",flopsNames[i],countFlops[i]);
    c+=countFlops[i];
  }
  fprintf(fp,"  Total   %9"PRId64"\n",c);

  long calls=0,execs=0, time=0;
  fprintf(fp,"Group   calls    execs   time [us]\n");
  //      <234> <2345678><2345678><23456789>
  for (int i=0;i<nGroups;i++) {
    calls+=countCallGroup[i];
    execs+=countExecuteGroup[i];
    time+=timeExecuteGroup[i];
    fprintf(fp,"%5d %9ld%9ld%14.0f\n",
	    i,countCallGroup[i],countExecuteGroup[i],timeExecuteGroup[i]*1e6/CLOCKS_PER_SEC);
  }
  fprintf(fp,"Total %9ld%9ld%14.0f\n",calls,execs,time*1e6/CLOCKS_PER_SEC);

  calls=0;time=0;
  fprintf(fp,"Set                                                                calls      time [us]\n");
  //        <234567890123456789012345678901234567890123456789> <2345678><23456789>
  for (int i=0;i<nSets;i++) {
    calls+=countCallSet[i];
    time+=timeExecuteSet[i];
    fprintf(fp,"  %-60s %9ld%14.0f\n",setNames[i],countCallSet[i],timeExecuteSet[i]*1e6/CLOCKS_PER_SEC);
  }
  fprintf(fp,"  Total                                                        %9ld%14.0f\n",
	 calls,time*1e6/CLOCKS_PER_SEC);

  calls=0;time=0;
  fprintf(fp,"Get                                                                calls      time [us]\n");
  //        <234567890123456789012345678901234567890123456789> <2345678><23456789>
  for (int i=0;i<nGets;i++) {
    calls+=countCallGet[i];
    time+=timeExecuteGet[i];
    fprintf(fp,"  %-60s %9ld%14.0f\n",getNames[i],countCallGet[i],timeExecuteGet[i]*1e6/CLOCKS_PER_SEC);
  }
  fprintf(fp,"  Total                                                        %9ld%14.0f\n",
	  calls,time*1e6/CLOCKS_PER_SEC);

  calls=0;time=0;
  fprintf(fp,"Copy                                                               calls      time [us]\n");
  //        <234567890123456789012345678901234567890123456789> <2345678><23456789>
  for (int i=0;i<nCopies;i++) {
    calls+=countCallCopy[i];
    time+=timeExecuteCopy[i];
    fprintf(fp,"  %-60s %9ld%14.0f\n",copyNames[i],countCallCopy[i],timeExecuteCopy[i]*1e6/CLOCKS_PER_SEC);
  }
  fprintf(fp,"  Total                                                        %9ld%14.0f\n",
	  calls,time*1e6/CLOCKS_PER_SEC);
}

EXPORT void profilingView(char * filename)
{
  FILE *fp=fopen(filename,"w");
  profilingViewFP(fp);
  fclose(fp);
  profilingViewFP(stdout);
}

EXPORT void profilingView0()
{
  FILE *fp=fopen("profileView.profile","w");
  profilingViewFP(fp);
  fclose(fp);
}
