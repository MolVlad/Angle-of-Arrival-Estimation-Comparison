# Angle-of-Arrival-Estimation-Comparison

Code for investigation of different smoothing modifications of the MUSIC algorithm for angle of arrival estimation. The link to the whole dataset can found
[here](http://wireless.iitp.ru/dataset-for-wi-fi-localization/). Dataset is gathered via USRP-based [testbed](https://github.com/MolVlad/LabVIEW-802.11-Angle-of-Arrival-Estimation).

## Files

- data/csi.mat -- channel state information from one experiment ([dataset](http://wireless.iitp.ru/dataset-for-wi-fi-localization/))
- data/samples.mat -- time domain IQ-samples from one experiment ([dataset](http://wireless.iitp.ru/dataset-for-wi-fi-localization/))
- data/snr.mat -- information about SNR from one experiment ([dataset](http://wireless.iitp.ru/dataset-for-wi-fi-localization/))
- calibration.mat -- calibration data ([info](http://wireless.iitp.ru/dataset-for-wi-fi-localization/))
- run_experiment.m -- main file for running experiments
- define_parameters.m -- file with settings algorithm parameters such as frequency, number of samples, smoothing type etc.
- tm.m -- implementation of time domain MUSIC with smoothing modifications
- fm.m -- implementation of frequency domain MUSIC with smoothing modifications
- findstable.m -- auxiliary function for finding the most stable peak on the graph

In case of any use, please cite original paper:

@article{molodtsov2021experimental,
  title={Experimental Study of Smoothing Modifications of the MUSIC Algorithm for Direction of Arrival Estimation in Indoor Environments},
  author={Molodtsov, Vladislav and Kureev, Aleksey and Khorov, Evgeny},
  journal={IEEE Access},
  volume={9},
  pages={153767--153774},
  year={2021},
  publisher={IEEE}
}
