#include "RecHitCalibration.hh"
#include <iostream>


RecHitCalibration :: RecHitCalibration(){
    float weights[53] =  {
                       0.0,   // there is no layer zero
                       8.603,  // Mev
                       8.0675,
                       8.0675,
                       8.0675,
                       8.0675,
                       8.0675,
                       8.0675,
                       8.0675,
                       8.0675,
                       8.9515,
                       10.135,
                       10.135,
                       10.135,
                       10.135,
                       10.135,
                       10.135,
                       10.135,
                       10.135,
                       10.135,
                       11.682,
                       13.654,
                       13.654,
                       13.654,
                       13.654,
                       13.654,
                       13.654,
                       13.654,
                       38.2005,
                       55.0265,
                       49.871,
                       49.871,
                       49.871,
                       49.871,
                       49.871,
                       49.871,
                       49.871,
                       49.871,
                       49.871,
                       49.871,
                       62.005,
                       83.1675,
                       92.196,
                       92.196,
                       92.196,
                       92.196,
                       92.196,
                       92.196,
                       92.196,
                       92.196,
                       92.196,
                       92.196,
                       46.098
    };
    
    float correction[3] = {1.132,1.092,1.084};
    float noises[3]     = {2100.0, 2100.0, 1600.0};  // 100,200,300 um (in electrons)
    float fcpermip[3]   = {1.25, 2.57, 3.88};  // 100um, 200um, 300um

    for (unsigned i=0; i<53; i++) dEdX_weights[i] = weights[i];
    for (unsigned i=0; i<3; i++) thicknessCorrection[i] = correction[i];
    for (unsigned i=0; i<3; i++) nonAgedNoises[i] = noises[i];
    for (unsigned i=0; i<3; i++) fCPerMIP[i] = fcpermip[i];
    
    fC_per_ele = 1.6020506e-4;
    noise_MIP = 1.0/7.0; //expectation based on latest SiPM performance
    }
    

float RecHitCalibration :: MeVperMIP(int layer, int thicknessIndex){
    if (layer > 40){
        return dEdX_weights[layer];
    }
    else{
        return dEdX_weights[layer]/thicknessCorrection[thicknessIndex];
    }
}

float RecHitCalibration :: MIPperGeV(int layer, int thicknessIndex){
    return 1000./MeVperMIP(layer, thicknessIndex);
}

float RecHitCalibration :: sigmaNoiseMIP(int layer, int thicknessIndex){
    if (layer > 40){
        // for BH, sigmaNoiseMIP = noise_MIP
        return noise_MIP;
    }
    else{
        return fC_per_ele *nonAgedNoises[thicknessIndex] / fCPerMIP[thicknessIndex];
    }
}

float RecHitCalibration :: sigmaNoiseMeV(int layer, int thicknessIndex){
    return sigmaNoiseMIP(layer, thicknessIndex) * MeVperMIP(layer, thicknessIndex);
}
