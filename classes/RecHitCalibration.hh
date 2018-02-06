#ifndef RECHITCALIBRATION_H
#define RECHITCALIBRATION_H

class RecHitCalibration{

    private:
        float dEdX_weights[53];
        float thicknessCorrection[3];
        float nonAgedNoises[3];  // 100,200,300 um (in electrons)

        // https://github.com/cms-sw/cmssw/blob/CMSSW_9_3_X/RecoLocalCalo/HGCalRecProducers/python/HGCalUncalibRecHit_cfi.py//L25
        float fCPerMIP[3];  // 100um, 200um, 300um

        // https://github.com/cms-sw/cmssw/blob/CMSSW_9_3_X/SimCalorimetry/HGCalSimProducers/python/hgcalDigitizer_cfi.py//L127
        float noise_MIP; //expectation based on latest SiPM performance
        float fC_per_ele;

    public:
        
	// constructor
	RecHitCalibration();
	float MeVperMIP(int layer, int thicknessIndex);
	float MIPperGeV(int layer, int thicknessIndex);
	float sigmaNoiseMIP(int layer, int thicknessIndex);
	float sigmaNoiseMeV(int layer, int thicknessIndex);

};

#endif 
