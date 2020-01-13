/*
 ============================================================================
 Name        : basisNquad.c
 Author      : Yukiko Shimizu
 Version     : 1.0
 Copyright   : Don't Touch. Hehe
 Description : Error Estimation for the Lorenz System
 ============================================================================
 */
#include "basisfunc.h"

double phi(int r, int k, double xi){  
  int i;  //initialization for iteration purposes
  double phi_k_xi= 1.0;   //Initialize to 1 for producting
  double xi_[r+1];    //Saves the xi value for each phi node
  //----------------------------------------------------------------------------
  // Implementation 
  //----------------------------------------------------------------------------
  for (i=0; i<r+1; i++)  //Calculates up to r+1 for each dt window
    xi_[i]=xiL(r, i);  //Calculates the all the xi nodes
  for (i=0; i<r+1; i++){  //Calculates Lagrange Basis function
    if (r>0){
      if (i!=k)  //When the phi of interest does not equal the node number
        phi_k_xi=phi_k_xi*(xi-xi_[i])/(xi_[k]-xi_[i]);
    }else{
      phi_k_xi = 1;
    }
  }
  return phi_k_xi;
}


double dphi(int r, int k, double xi){
  int i, j;
  double dphi_k_xi=0;   //Initialize to 1 for producting
  double xi_[r+1];    //Saves the xi value for each phi node
  double temp[r+1];  //Array that holds intermediate calculations for phis
  double temp2[r+1]; //Array that holds intermediate calculations for dphis
  //===========================Implementation==================================
  //Calculates phi, dphi, and nodes(xiL) in reference space
  for (i=0; i<r+1; i++)  //Calculates up to r+1 for each dt window
    xi_[i]=xiL(r, i);  //Calculates the all the xi nodes
  if (r>0){
    for (i=0; i<r+1; i++){ //Calculates parts of Lagrange Basis function
      if (i!=k){      //When the phi of interest does not equal the node number
        //Saves parts of the calculation to temp for further calculation
        temp[i]=(xi-xi_[i])/(xi_[k]-xi_[i]);
      }else{
        temp[i]=0;  //Not used when k==j since not part of the Lagrange Basis
      }
    }
    //Calculation of the intermediate values for dphi
    for (i=0; i<r+1; i++){
      if (i!=k)  //when phi number does not equal node number
        //Derivative of the Lagrange Basis function parts
        temp2[i]=1.0/(xi_[k]-xi_[i]); //temp2 saves intermediate values
    }
  }
  //Calculation of final dphis with a final summation at the end
  for (i=0; i<r+1; i++){
    if (r>0){
      if (i!=k){  //when phi numebr does not equal node number
        for (j=0; j<r+1; j++){ //For loop to keep track of the chain rule
          if (j!=i && j!=k && r!=1 )
            temp2[i]=temp2[i]*temp[j];  //First performs chain rule for deriva 
        }
        //Now perform last part of chain rule to get the final dphi values
        dphi_k_xi+=temp2[i];
      }
    }else{
      dphi_k_xi=0;
    }
  }
  return dphi_k_xi;
}


double d2phi(int p, int k, double xi){
  int i, j, n;    //initialization for iteration purposes
  double xi_[p+1];
  double temp[p+1];
  double dtemp[p+1];   //tmeporary values evaluated for derivative of basis
  double indtemp;  //initialization to 1 for product sum
  double d2phi_k_xi = 0;
  for (i=0; i<p+1; i++)
    xi_[i] = xiL(p,i);
  if (p>1){
    //Calculates basis function phi for all vlaue of i given k
    for (i=0; i<p+1; i++){
      if (i!=k){
        temp[i] = (xi-xi_[i])/(xi_[k]-xi_[i]);
      }else{
        temp[i] = 0;
      }
    }
    //Calculates the derivative of the basis functions dphi
    //for all value of i given k
    for (i=0; i<p+1; i++){
      if (i!=k)
        dtemp[i] = 1.0/(xi_[k]-xi_[i]);
    }
  }

  if (p>1){
    for (i=0; i<p+1; i++){  //iterate through p+1
      for (j=i; j<p+1; j++){
        indtemp =1;
        if (i!=k && j!=k && i!=j){
          for (n=0; n<p+1; n++){
            if (i!= n && j!= n && n!=k)
              indtemp *=temp[n];
          }
          indtemp *= dtemp[i]*dtemp[j];
          d2phi_k_xi +=2.0*indtemp;
        }
      }
    }
  }else{
    d2phi_k_xi = 0;
  }
  return d2phi_k_xi;
}          


double d3phi(int p, int k, double xi){
  int i, j, m, n;  //initialization for iteraiton purposes
  double xi_[p+1];
  double temp[p+1];
  double dtemp[p+1]; //temporary values evlatued for derivative of basis
  double indtemp;
  double d3phi_k_xi =0;

  for (i=0; i<p+1; i++)
    xi_[i] = xiL(p,i);
  if (p>1){
    //Calculates basis function phi for all vlaue of i given k
    for (i=0; i<p+1; i++){
      if (i!=k){
        temp[i] = (xi-xi_[i])/(xi_[k]-xi_[i]);
      }else{
        temp[i] = 0;
      }
    }
    //Calculates the derivative of the basis functions dphi
    //for all value of i given k
    for (i=0; i<p+1; i++){
      if (i!=k)
        dtemp[i] = 1.0/(xi_[k]-xi_[i]);
    }
  }

  if (p>2){
    for (i=0; i<p+1; i++){
      for (j=i; j<p+1; j++){
        for (m=j; m<p+1; m++){
          indtemp = 1;
          if (i!=k && j!=k && m!=k && i!=j && j!=m &&i!=m){
            for (n=0; n<p+1; n++){
              if (i!= n && j!= n && m!=n && n!=k)
                indtemp *=temp[n];
            }
            indtemp *= dtemp[i]*dtemp[j]*dtemp[m];
            d3phi_k_xi +=6.0*indtemp;
          }
        }
      }
    }
  }else{
    d3phi_k_xi = 0;
  }
  return d3phi_k_xi;
}                  

double d4phi(int p, int k, double xi){
  int i, j, m, n; //initialization for iteration purposes
  int q; //initialixation for iteration purposes

  double xi_[p+1];
  double temp[p+1];
  double dtemp[p+1]; //temporary values evlatued for derivative of basis
  double indtemp;
  double d4phi_k_xi = 0;

  for (i=0; i<p+1; i++)
    xi_[i] = xiL(p,i);
  if (p>1){
    //Calculates basis function phi for all vlaue of i given k
    for (i=0; i<p+1; i++){
      if (i!=k){
        temp[i] = (xi-xi_[i])/(xi_[k]-xi_[i]);
      }else{
        temp[i] = 0;
      }
    }
    //Calculates the derivative of the basis functions dphi
    //for all value of i given k
    for (i=0; i<p+1; i++){
      if (i!=k)
        dtemp[i] = 1.0/(xi_[k]-xi_[i]);
    }
  }
  if (p>3){
    for (i=0; i<p+1; i++){
      for (j=i; j<p+1; j++){
        for (m=j; m<p+1; m++){
          for(q=m; q<p+1; q++){
            indtemp = 1;
            if (i!=k && j!=k && m!=k && q!=k
                && i!=j && j!=m &&i!=m && i!=q && j!=q && m!= q){
              for (n=0; n<p+1; n++){
                if (i!= n && j!= n && m!=n && q!=n && n!=k)
                  indtemp *=temp[n];
              }
              indtemp *= dtemp[i]*dtemp[j]*dtemp[m]*dtemp[q];
              d4phi_k_xi +=24.0*indtemp;
            }
          }
        }
      }
    }
  }else{
    d4phi_k_xi = 0;
  }
  return d4phi_k_xi;
}                
