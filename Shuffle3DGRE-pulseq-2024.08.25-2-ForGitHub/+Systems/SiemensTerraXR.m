function [SysNormal, SystemsStructArray, SysNormalIdx] = SiemensTerraXR()

% [SysNormal, SysStructArray, SysNormalIdx] = SiemensTerraXR()
%
% -------------------------------------------------------------------------
% Terra (7T) system defintions
% -------------------------------------------------------------------------
%
% Returns system definitions, to be used by pulseq, for the Siemens Terra.
% (Also know as 'Terra-XR', including 'Terra-XR_MNO', 'Terra-XR_8Tx', and
% 'Terra-XR_16Tx'. Maybe others as well in the future.)
%
% Outputs:
% --------
%
% - SysNormal - The default system to be used. Includes "standard"
%   gradient amplitude and slew rate limits.
% - SysStructArray - An array of all(!) systems available (including the
%   "normal" one). Each element SysStructArray(i) includes, at least two
%   fields:
%   > SysStructArray(i).System - a pulseq appropriate system definition.
%   > SysStructArray(i).Name - a name/description of the system.
%   The systems in the array are expected to be ordered from "safest"/least
%   demanding to the "worst"/most demanding system.
% - SysNormalIdx - The index of SysNormal within SysStructArray ;



%% Base system (without limits on gradient amplitude/rise)

% Since several gradient limits (fast/normal/whisper) are possible, we
% first define limits without gradients and then base gradient specific
% limits on these.
%
% NOTE: Text in parantheses, e.g., "NominalB0", is the description in the
%       output of the 'imprint' command (in IDEA or on the scanner), the
%       value of which was used to set the parameters given below.


SystemNoGrad.B0 = 6.980936 ; % [T] - "NominalB0"
SystemNoGrad.gamma = 42575575 ; % [Hz/T] - "Larmor Constant 1H"

SystemNoGrad.rfDeadTime = 100e-6 ; % [s] - "RFCI CoilCtrlLead time"
SystemNoGrad.rfRingdownTime = 20e-6 ; % [s] - "RFCI CoilCtrlHold time"
% Note: actual time between RF and ADC, ADC and RF, or ADC and ADC depends
%       on the dwell time of the first "action" of the two that is run.
%       In IDEA this can be found by running one of the following: 
%         getMinDurationBetweenFrequencyPhaseEvents()
%         getMinDurationBetweenReadoutAndRFPuls()
%         getMinDurationBetweenReadoutAndReadout()
%         getMinDurationBetweenFrequencyPhaseEvents()
%      If the two concuring ADCs have the same dwell time is the same, the
%      required delay is zero.
%      We will use a value of 10us so that a block with only an ADC in it
%      (noise adjust scan) will not fail. If we have two consective ADCs
%      without a delay between them, this might be a problem, unless pulseq
%      solves this somehow (e.g. concatenates ADCs when the requested ADC
%      is too long).
SystemNoGrad.adcDeadTime = 10e-6 ; % [s] see note above(!)

SystemNoGrad.adcRasterTime = 100e-9 ; % [s]
SystemNoGrad.rfRasterTime = 1e-6 ; % [s]
SystemNoGrad.gradRasterTime = 10e-6 ; % [s]
SystemNoGrad.blockDurationRaster = 10e-6 ; % [s]

%% Grdaient "systems"/modes

SysCounter = 0 ;


% System Whisper
% --------------

SysCounter = SysCounter + 1 ;

% Short Name, for convenience
System = SystemNoGrad ;
System.maxGrad = mr.convert(22,'mT/m', ... % "GPA GradMaxAmpl" (#3)
                           'Hz/m', ...
                            'gamma', SystemNoGrad.gamma) ; % [Hz/m] 
System.maxSlew = mr.convert(100,'T/m/s', ... % "GPA GradMaxSlewRate" (#3)
                            'Hz/m/s', ...
                            'gamma', SystemNoGrad.gamma) ; % [Hz/m/s] 
System.riseTime = [] ; % [s] Not really used anywhere.

% update SystemsStructArray
SystemsStructArray(SysCounter).System = System ;
SystemsStructArray(SysCounter).Name = 'Whisper' ;


                                                   


% System Normal
% -------------

SysCounter = SysCounter + 1 ;

% Short Name, for convenience
System = SystemNoGrad ;
System.maxGrad = mr.convert(40,'mT/m', ... % "GPA GradMaxAmpl" (#2)
                           'Hz/m', ...
                            'gamma', SystemNoGrad.gamma) ; % [Hz/m] 
System.maxSlew = mr.convert(200,'T/m/s', ... % "GPA GradMaxSlewRate" (#2)
                            'Hz/m/s', ...
                            'gamma', SystemNoGrad.gamma) ; % [Hz/m/s] 
System.riseTime = [] ; % [s] Not really used anywhere.

% update SystemsStructArray
SystemsStructArray(SysCounter).System = System ;
SystemsStructArray(SysCounter).Name = 'Normal' ;

% Set as "normal" system:
SysNormal = System ;
SysNormalIdx = SysCounter ;



% SystemFast
% ----------

SysCounter = SysCounter + 1 ;

% Short Name, for convenience
System = SystemNoGrad ;
System.maxGrad = mr.convert(42,'mT/m', ... % "GPA GradMaxAmpl" (#1)
                           'Hz/m', ...
                            'gamma', SystemNoGrad.gamma) ; % [Hz/m] 
System.maxSlew = mr.convert(200,'T/m/s', ... % "GPA GradMaxSlewRate" (#1)
                            'Hz/m/s', ...
                            'gamma', SystemNoGrad.gamma) ; % [Hz/m/s] 
System.riseTime = [] ; % [s] Not really used anywhere.

% update SystemsStructArray
SystemsStructArray(SysCounter).System = System ;
SystemsStructArray(SysCounter).Name = 'Fast' ;




end