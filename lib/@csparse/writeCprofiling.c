/*
  Copyright 2012-2017 Joao Hespanha

  This file is part of Tencalc.

  TensCalc is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  TensCalc is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with TensCalc.  If not, see <http://www.gnu.org/licenses/>.
*/

void profilingView(char * filename)
{
  int nGroups=sizeof(countCallGroup)/sizeof(countCallGroup[0]);
  int nGets  =sizeof(getNames)/sizeof(getNames[0]);
  int nSets  =sizeof(setNames)/sizeof(setNames[0]);
  int nCopies=sizeof(copyNames)/sizeof(copyNames[0]);
  int nFlops =sizeof(flopsNames)/sizeof(flopsNames[0]);

  FILE *fp=fopen(filename,"w");
  
  printf("Operation    Count\n");
  fprintf(fp,"Operation    Count\n");
  int64_t c=0;
  for (int i=0;i<nFlops;i++) {
    printf("%-10s%9"PRId64"\n",flopsNames[i],countFlops[i]);
    fprintf(fp,"%-10s%9"PRId64"\n",flopsNames[i],countFlops[i]);
    c+=countFlops[i];
  }
  printf("  Total   %9"PRId64"\n",c);
  fprintf(fp,"  Total   %9"PRId64"\n",c);
  
  long calls=0,execs=0, time=0;
  printf("Group   calls    execs   time [us]\n");
  //      <234> <2345678><2345678><23456789>
  for (int i=0;i<nGroups;i++) {
    calls+=countCallGroup[i];
    execs+=countExecuteGroup[i];
    time+=timeExecuteGroup[i];
    printf("%5d %9ld%9ld%14.0f\n",
	   i,countCallGroup[i],countExecuteGroup[i],timeExecuteGroup[i]*1e6/CLOCKS_PER_SEC);
    fprintf(fp,"%5d %9ld%9ld%14.0f\n",
	    i,countCallGroup[i],countExecuteGroup[i],timeExecuteGroup[i]*1e6/CLOCKS_PER_SEC);
  }
  printf("Total %9ld%9ld%14.0f\n",calls,execs,time*1e6/CLOCKS_PER_SEC); 
  fprintf(fp,"Total %9ld%9ld%14.0f\n",calls,execs,time*1e6/CLOCKS_PER_SEC); 

  calls=0;time=0;
  printf("Set                                                                calls      time [us]\n");
  fprintf(fp,"Set                                                                calls      time [us]\n");
  //        <234567890123456789012345678901234567890123456789> <2345678><23456789>
  for (int i=0;i<nSets;i++) {
    calls+=countCallSet[i];
    time+=timeExecuteSet[i];
    printf("  %-60s %9ld%14.0f\n",setNames[i],countCallSet[i],timeExecuteSet[i]*1e6/CLOCKS_PER_SEC);
    fprintf(fp,"  %-60s %9ld%14.0f\n",setNames[i],countCallSet[i],timeExecuteSet[i]*1e6/CLOCKS_PER_SEC);
  }
  printf("  Total                                                        %9ld%14.0f\n",
	 calls,time*1e6/CLOCKS_PER_SEC); 
  fprintf(fp,"  Total                                                        %9ld%14.0f\n",
	 calls,time*1e6/CLOCKS_PER_SEC); 

  calls=0;time=0;
  printf("Get                                                                calls      time [us]\n");
  fprintf(fp,"Get                                                                calls      time [us]\n");
  //        <234567890123456789012345678901234567890123456789> <2345678><23456789>
  for (int i=0;i<nGets;i++) {
    calls+=countCallGet[i];
    time+=timeExecuteGet[i];
    printf("  %-60s %9ld%14.0f\n",getNames[i],countCallGet[i],timeExecuteGet[i]*1e6/CLOCKS_PER_SEC);
    fprintf(fp,"  %-60s %9ld%14.0f\n",getNames[i],countCallGet[i],timeExecuteGet[i]*1e6/CLOCKS_PER_SEC);
  }
  printf("  Total                                                        %9ld%14.0f\n",
	 calls,time*1e6/CLOCKS_PER_SEC); 
  fprintf(fp,"  Total                                                        %9ld%14.0f\n",
	  calls,time*1e6/CLOCKS_PER_SEC); 

  calls=0;time=0;
  printf("Copy                                                               calls      time [us]\n");
  fprintf(fp,"Copy                                                               calls      time [us]\n");
  //        <234567890123456789012345678901234567890123456789> <2345678><23456789>
  for (int i=0;i<nCopies;i++) {
    calls+=countCallCopy[i];
    time+=timeExecuteCopy[i]; 
    printf("  %-60s %9ld%14.0f\n",copyNames[i],countCallCopy[i],timeExecuteCopy[i]*1e6/CLOCKS_PER_SEC);
    fprintf(fp,"  %-60s %9ld%14.0f\n",copyNames[i],countCallCopy[i],timeExecuteCopy[i]*1e6/CLOCKS_PER_SEC);
  }
  printf("  Total                                                        %9ld%14.0f\n",
	 calls,time*1e6/CLOCKS_PER_SEC); 
  fprintf(fp,"  Total                                                        %9ld%14.0f\n",
	  calls,time*1e6/CLOCKS_PER_SEC);
  fclose(fp);
}
