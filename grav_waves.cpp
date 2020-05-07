// Matt Walsh

/* I found four matches of the wave template in the hanford data and five in the livingston data but only one overlapped, meaning there was only one gravitational wave. The overlap was recorded at about t=0.9514 seconds. As i said above, both data sets had multiple points where the chi2 value of the template vs the data was very low, howerver only one point had a low chi2 value at the same time point in both data sets */

/*My approach to this problem and explanation of code: I first read all the files and determined how many entries each had and filled the histogram to find the mean to use for sigma in my chi2 method. I then created arrays to put all the data into so that i could compare them easily, and I used the number of entries i found from before to use as the size of the array. I reread each file to fill each of the arrays. I then used for loops to find a chi2 model of the data vs the templates by moving the template along the data and comparing each point. I then graphed these chi2 values to find where the chi2 value dipped very low, meaning the data matched the template. I did this for both data sets and compared the two graphs to eachother to find where both data sets showed dips in their respective chi2 values*/

#include <stdio.h>
#include <stdlib.h>
#include "hist.hpp"

int main() {
  h1 hist;
  h1init(&hist, 100, -1, 1, "Strain Histogram");
  h1labels(&hist, "Strain * 1.e21", "Frequency");

  FILE* h_dat;
  FILE* l_dat;
  FILE* h_wave;
  FILE* l_wave;
  
  h_dat = fopen("hanford_waveform_complete.dat", "r");
  l_dat = fopen("livingston_waveform_complete.dat", "r");
  h_wave = fopen("hanford_GRmodel.dat", "r");
  l_wave = fopen("livingston_GRmodel.dat", "r");
  

  if (h_dat==NULL||l_dat==NULL){
    printf("Error: file not found!\n");
    return -1;
  }

  int status=0;
  double time, l_time, h_time, strain;
  int num_h = 0; // number of entries in h file
  int num_l = 0; // number of entries in l file
  int num_w_h = 0; // number of entries in wave template h
  int num_w_l = 0; // number of entries in wave template l

  char titles[1000]; //these lines will skip past the titles in the data files
  fgets(titles, 1000, h_dat); //titles must be less than 1000 chars long
  fgets(titles, 1000, l_dat);
  fgets(titles, 1000, h_wave);
  fgets(titles, 1000, l_wave);

  //read though files to fill histogram and count entries
  while (status!=EOF) {
    status = fscanf(h_dat, "%lf %lf", &time, &strain);

    if (status == 2) {
      h1fill(&hist, strain);
      num_h++;
    }
  }

  status = 0;

  while (status!=EOF) {
    status = fscanf(l_dat, "%lf %lf", &time, &strain);

    if (status == 2) {
      h1fill(&hist, strain);
      num_l++;
    }
  }

  status = 0;

  while (status!=EOF) { //to find total time of h wave template
    status = fscanf(h_wave, "%lf %lf", &h_time, &strain);
    if (status==2) {
      num_w_h++;
    }
  }
  
  status = 0;

  while (status!=EOF) { //to find total time of l wave template
    status = fscanf(l_wave, "%lf %lf", &l_time, &strain); 
    if (status==2) {
      num_w_l++;
    }
  }

  //reset where we are reading in file
  rewind(h_dat);
  rewind(l_dat);
  rewind(h_wave);
  rewind(l_dat);
  
  //pull data from hist
  int entries;
  double mean;
  double std_dv;
  h1stats(&hist, &entries, &mean, &std_dv);

  double l_arr[2][num_l]; //array to hold l data
  double h_arr[2][num_h]; // array to hold h data
  double l_warr[2][num_w_l]; // array to hold l wave template data
  double h_warr[2][num_w_h]; // array to hold h wave template data

  h_dat = fopen("hanford_waveform_complete.dat", "r");
  l_dat = fopen("livingston_waveform_complete.dat", "r");
  h_wave = fopen("hanford_GRmodel.dat", "r");
  l_wave = fopen("livingston_GRmodel.dat", "r");

  fgets(titles, 1000, h_dat); //titles must be less than 1000 chars long
  fgets(titles, 1000, l_dat);
  fgets(titles, 1000, h_wave);
  fgets(titles, 1000, l_wave);

  double th, tl, twh, twl; //time values for all data sets


  //fill arrays
  status = 0;
  int count = 0;
  while (status!=EOF) {
    status = fscanf(h_dat, "%lf %lf", &th, &strain);
    if (status == 2) {
      h_arr[0][count]=th;
      h_arr[1][count]=strain;
      count++;
    }
  }

  status = 0;
  count = 0;

  while (status!=EOF) {
    status = fscanf(l_dat, "%lf %lf", &tl, &strain);
    if (status == 2) {
      l_arr[0][count]=tl;
      l_arr[1][count]=strain;
      count++;
    }
  }

  status = 0;
  count = 0;
  
while (status!=EOF) {
    status = fscanf(h_wave, "%lf %lf", &twh, &strain);
    if (status == 2) {
      h_warr[0][count]=twh;
      h_warr[1][count]=strain;
      count++;
    }
  }

  status = 0;
  count = 0;

while (status!=EOF) {
    status = fscanf(l_wave, "%lf %lf", &twl, &strain);
    if (status == 2) {
      l_warr[0][count]=twl;
      l_warr[1][count]=strain;
      count++;
    }
  }
// arrays are now filled

 FILE* hout = fopen("ligo_out_h.dat", "w");
  double chi2=0;
  for (int a = 0; a < num_h-num_w_h; a++) {
    //first loop goes thorugh complete waveform
    chi2 = 0;
    for (int b = 0; b < num_w_h; b++) {
      //second loop compares wave template so section of complete waveform
      //calculate chi squared
      chi2 += ((h_arr[1][a+b]-h_warr[1][b])*(h_arr[1][a+b]-h_warr[1][b]))/(mean*mean);
    }
    chi2/=10000000;//not sure why i need this but chi2 values were 2e7 ish
    fprintf(hout, "%lf %lf\n", h_arr[0][a], chi2);
  }

  FILE* lout = fopen("ligo_out_l.dat", "w");
  for (int a = 0; a < num_l-num_w_l; a++) {
    //first loop goes thorugh complete waveform
    chi2 = 0;
    for (int b = 0; b < num_w_l; b++) {
      //second loop compares wave template so section of complete waveform
      //calculate chi squared
      chi2 += ((l_arr[1][a+b]-l_warr[1][b])*(l_arr[1][a+b]-l_warr[1][b]))/(mean*mean);
    }
    chi2/=10000000;//not sure why i need this but chi2 values were 2e7 ish
    fprintf(lout, "%lf %lf\n", l_arr[0][a], chi2);
  }
  
  
  
}
