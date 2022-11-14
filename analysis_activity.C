#include "TH1.h"
#include "TF1.h"
#include "TGraph2DErrors.h"
#include "TMath.h"


 
void analysis()
{

  ifstream inputfile;
  inputfile.open("gplevh.dat", ifstream::in);
  if( !inputfile ) {
    cerr << "Error opening input stream" << endl;
    return;
  }  


  double dx = 0.1;
  double dy = 0.1;
  double norm = 27.027027027;
   const int np = 441;
   double x[np];
   double ex[np];
   double y[np];
   double ey[np];
   double act[np];
   double eact[np];

   double p0, p1, p2, p3;
   int i=0;
   
   while( !inputfile.eof() ){
 
     inputfile >> p0 >> p1 >> p2 >> p3;
     if (p2>0) {
       cout << p0 << " " << p1 << " " << p2 << " " << p3 << endl;
       x[i] = p0 +dx/2;
       y[i] = p1 +dy/2;
       ex[i] = dx/2;
       ey[i] = dy/2;
       act[i]    = p2 * norm ;
       eact[i]   = p3 *act[i] / 100;
       i=i+1;
     }
   }
   inputfile.close();

   TGraph2DErrors *gra = new TGraph2DErrors(i,x,y,act,ex,ey,eact);
   gra->SetTitle("Chromium 51 activity map 1month@10#muA;X(cm);Y(cm);pCi/cm^{3}");
   gra->GetXaxis()->SetTitleOffset(1.5);
   gra->GetYaxis()->SetTitleOffset(2.5);
   gra->GetZaxis()->SetTitleOffset(2.5);
   //   gra->SetMinimum(1.0E5);

   TCanvas * c1 = new TCanvas();
   c1->SetLogz();
   gra->Draw("TRI2ZERR");
   c1->Print("graph_test.pdf");
   

}


void analysis_sym()
{

  ifstream inputfile;
  inputfile.open("gplevh.dat", ifstream::in);
  if( !inputfile ) {
    cerr << "Error opening input stream" << endl;
    return;
  }  


  double dx = 0.1;
  double dy = 0.1;
  double norm = 27.027027027;
  const int np = 441;
  double x[np];
  double ex[np];
  double y[np];
  double ey[np];
  double xsym[np];
  double exsym[np];
  double ysym[np];
  double eysym[np];
  double act[np];
  double eact[np];
  double actsym[np];
  double eactsym[np];

  double p0, p1, p2, p3;
  int i=0;
   
  while( !inputfile.eof() ){
 
    inputfile >> p0 >> p1 >> p2 >> p3;
    if (p2>0) {
      cout << p0 << " " << p1 << " " << p2 << " " << p3 << endl;
      x[i] = TMath::Abs(p0+dx/2);
      y[i] = TMath::Abs(p1+dy/2);
      ex[i] = dx/2;
      ey[i] = dy/2;
      act[i]    = p2 * norm ;
      eact[i]   = p3 *act[i] / 100;
      i=i+1;
    }
  }
  inputfile.close();
  int j=0;
  int k=0;
  int ihave =0;
  for (j=0; j<i; j++) {
    ihave = 0;
    cout << "At " << j << " of total of " << i << endl;
    for (int m=0; m<k ; m++) {
      if (x[j]==xsym[m] && y[j]==ysym[m] && eact[j]> 0.0) { // I have already this cohordinate
	actsym[m] = (actsym[m]/pow(eactsym[m],2) + act[j]/pow(eact[j],2))/(1./pow(eactsym[m],2) + 1./pow(eact[j],2));
	eactsym[m] = pow(1./(1./pow(eactsym[m],2) + 1./pow(eact[j],2)),0.5); // weighted average and error
	ihave = 1;
	cout << "Double: x=" << x[j] << " y=" << y[j] << endl;
      }
      else if (y[j]==xsym[m] && x[j]==ysym[m] && eact[j]> 0.0 && x[j] != y[j] ) { // I have already this cohordinate using symmetry in x and y, but diagonal cannot count twice
	actsym[m] = (actsym[m]/pow(eactsym[m],2) + act[j]/pow(eact[j],2))/(1./pow(eactsym[m],2) + 1./pow(eact[j],2));
	eactsym[m] = pow(1./(1./pow(eactsym[m],2) + 1./pow(eact[j],2)),0.5); // weighted average and error
	ihave = 1;
	cout << "Double: x=" << x[j] << " y=" << y[j] << endl;
      }
    }
    if (ihave == 0 && eact[j] > 0) {// I don't have the point and it is > 0
      xsym[k] = x[j];
      ysym[k] = y[j];
      exsym[k] = ex[j];
      eysym[k] = ey[j];
      actsym[k] = act[j];
      eactsym[k] = eact[j];
      k++;
    }
  }

  double dummy;
  for (j=0; j<k; j++) {
    if (ysym[j] > xsym[j]) {
      xsym[j] = dummy;
      xsym[j] = ysym[j];
      ysym[j] = dummy;
    }
  }

  for (j=0; j<k; j++) {
    xsym[j+k] = ysym[j];
    ysym[j+k] = xsym[j];
    exsym[j+k] = exsym[j];
    eysym[j+k] = eysym[j];
    actsym[j+k] = actsym[j];
    eactsym[j+k] = eactsym[j];
    xsym[j+2*k] = -xsym[j];
    ysym[j+2*k] = ysym[j];
    exsym[j+2*k] = exsym[j];
    eysym[j+2*k] = eysym[j];
    actsym[j+2*k] = actsym[j];
    eactsym[j+2*k] = eactsym[j];
    xsym[j+3*k] = -ysym[j];
    ysym[j+3*k] = xsym[j];
    exsym[j+3*k] = exsym[j];
    eysym[j+3*k] = eysym[j];
    actsym[j+3*k] = actsym[j];
    eactsym[j+3*k] = eactsym[j];
    xsym[j+4*k] = xsym[j];
    ysym[j+4*k] = -ysym[j];
    exsym[j+4*k] = exsym[j];
    eysym[j+4*k] = eysym[j];
    actsym[j+4*k] = actsym[j];
    eactsym[j+4*k] = eactsym[j];
    xsym[j+5*k] = ysym[j];
    ysym[j+5*k] = -xsym[j];
    exsym[j+5*k] = exsym[j];
    eysym[j+5*k] = eysym[j];
    actsym[j+5*k] = actsym[j];
    eactsym[j+5*k] = eactsym[j];
    xsym[j+6*k] = -xsym[j];
    ysym[j+6*k] = -ysym[j];
    exsym[j+6*k] = exsym[j];
    eysym[j+6*k] = eysym[j];
    actsym[j+6*k] = actsym[j];
    eactsym[j+6*k] = eactsym[j];
    xsym[j+7*k] = -ysym[j];
    ysym[j+7*k] = -xsym[j];
    exsym[j+7*k] = exsym[j];
    eysym[j+7*k] = eysym[j];
    actsym[j+7*k] = actsym[j];
    eactsym[j+7*k] = eactsym[j];

  }
  
  TGraph2DErrors *gra = new TGraph2DErrors(8*k,xsym,ysym,actsym,exsym,eysym,eactsym);
  gra->SetTitle("Chromium 51 activity map 1month@10#muA;X(cm);Y(cm);pCi/cm^{3}");
  gra->GetXaxis()->SetTitleOffset(2.5);
  gra->GetYaxis()->SetTitleOffset(1.5);
  gra->GetZaxis()->SetTitleOffset(1.5);
  gra->SetMinimum(1.0E5);
  gra->SetMaximum(5.0E10);

  TCanvas * c1 = new TCanvas();
  c1->SetLogz();
  gra->SetMarkerStyle(20);
  gra->Draw("PTRI2ERRZ");
  c1->Print("graph_testsym.pdf");
   

}
