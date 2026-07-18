@echo off
rm synthese_observations.csv
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_07_16_10_12-180-0-SRT-Krakow-3m.fits.csv" 0 0 1.14 0 -a 6.3
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_07_20_45_12--175-0-SRT-Krakow-3m.fits.csv" 0 0 -1.4 0
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_07_24_44_12--170-0-SRT-Krakow-3m.fits.csv" 0 0 -4 0 -a 3.5
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_07_28_29_12--165-0-SRT-Krakow-3m.fits.csv" 0 0 -4 0 -a 10.1
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_07_32_19_12--160-0-SRT-Krakow-3m.fits.csv" 0 0 -8 0 -a 10.2
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_07_36_01_12--155-0-SRT-Krakow-3m.fits.csv" 0 0 -14 0
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_07_39_36_12--150-0-SRT-Krakow-3m.fits.csv" 0 0 -14 0 -r 3 -l cl
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_07_43_11_12--145-0-SRT-Krakow-3m.fits.csv" 0 100 -15.5 0 -l cl
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_07_46_50_12--140-0-SRT-Krakow-3m.fits.csv" -10 100 -15 0 -l cl -a 7.3
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_07_50_39_12--135-0-SRT-Krakow-3m.fits.csv" -20 100 -21 0 -l cl -a 49.2
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_07_54_15_12--130-0-SRT-Krakow-3m.fits.csv" -10 100 -21 0 -l cl
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_14_11_20_44_14--125-0-SRT-Krakow-3m.fits.csv" -10 100 -22.5 0 -l cl
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_14_11_16_59_14--120-0-SRT-Krakow-3m.fits.csv" -10 100 -24 0 -l cl -p 0.1
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_14_11_13_11_14--115-0-SRT-Krakow-3m.fits.csv" -10 100 -25 0 -l cl -a 0.7
REM python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_08_04_26_12-175-0-SRT-Krakow-3m.fits.csv"  -50 30 4 0 -l cl
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_14_12_22_10_14-175-0-SRT-Krakow-3m.fits.csv" -50 30 8 0
REM python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_08_08_00_12-170-0-SRT-Krakow-3m.fits.csv" -60 30 4 0 -l cr
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_14_12_18_22_14-170-0-SRT-Krakow-3m.fits.csv" -50 30 0 0
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_14_12_14_48_14-165-0-SRT-Krakow-3m.fits.csv" -80 30 -10 0 -a -17
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_14_12_11_07_14-160-0-SRT-Krakow-3m.fits.csv" -80 20 10 0
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_14_12_06_36_14-155-0-SRT-Krakow-3m.fits.csv" -90 40 10 0 -r 1,2 -a=-31.7,-18.6
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_14_12_02_52_14-150-0-SRT-Krakow-3m.fits.csv" -90 20 13 0
REM python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_08_14_15_12-145-0-SRT-Krakow-3m.fits.csv" -90 20 20 0 -p 0.2 -l cr -a -69.3
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_14_11_59_06_14-145-0-SRT-Krakow-3m.fits.csv" -90 20 19.5 0 -p 0.2 -l cr
REM python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_08_18_11_12-140-0-SRT-Krakow-3m.fits.csv" -100 10 18 0 -l cr -a=-72.4,-87.0,-19.3
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_14_11_30_22_14-140-0-SRT-Krakow-3m.fits.csv" -100 10 18 0 -l cr -a=-72.4,-87.0,-19.3
REM python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_08_21_43_12-135-0-SRT-Krakow-3m.fits.csv" -120 20 23 0 -l cr -p 0.5 -a -91.4
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_14_11_26_41_14-135-0-SRT-Krakow-3m.fits.csv" -120 20 26 0 -l cr
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_08_25_16_12-130-0-SRT-Krakow-3m.fits.csv" -120 20 22 0 -l cr -a -86.3
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_08_28_47_12-125-0-SRT-Krakow-3m.fits.csv" -100 20 20 0 -l cr -a -98
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_08_32_19_12-120-0-SRT-Krakow-3m.fits.csv" -80 20 23 0 -l cr -a=-96.0,-32.0
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_08_48_01_12-115-0-SRT-Krakow-3m.fits.csv" -110 20 22 0 -l cr -a=-100,-35.9
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_08_51_37_12-110-0-SRT-Krakow-3m.fits.csv" -100 20 25 0 -l cr -a=-96.4,-24.6
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_08_55_39_12-105-0-SRT-Krakow-3m.fits.csv" -100 20 28 0 -l cr -a -90.6
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_08_59_11_12-100-0-SRT-Krakow-3m.fits.csv" -100 30 28 0 -l cr
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_09_03_09_12-95-0-SRT-Krakow-3m.fits.csv" -100 30 32 0 -l cr -a=-44,-5.9
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_09_06_52_12-90-0-SRT-Krakow-3m.fits.csv" -100 30 25 0 -l cr -a -21.4
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_20_31_53_12-85-0-SRT-Krakow-3m.fits.csv" -100 40 27 0 -a -11.1
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_20_38_40_12-80-0-SRT-Krakow-3m.fits.csv" -100 40 25 0 -p 0.2
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_20_42_46_12-75-0-SRT-Krakow-3m.fits.csv" -100 40 25 0 -p 0.3 -a -27.9
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_20_46_28_12-70-0-SRT-Krakow-3m.fits.csv" -90 50 23 0 -p 0.2
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_20_51_05_12-65-0-SRT-Krakow-3m.fits.csv" -100 60 25 0 -a -12.7
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_20_54_54_12-60-0-SRT-Krakow-3m.fits.csv" -100 70 25 0 -p 0.35 -a 28.9
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_20_58_37_12-55-0-SRT-Krakow-3m.fits.csv" -90 70 19 0 -p 0.3
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_21_02_22_12-50-0-SRT-Krakow-3m.fits.csv" -80 90 20 0 -p 0.2
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_21_06_12_12-45-0-SRT-Krakow-3m.fits.csv" -80 100 20 0 -p 0.2
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_21_10_02_12-40-0-SRT-Krakow-3m.fits.csv" -80 100 20 0 -p 0.2
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_21_14_03_12-35-0-SRT-Krakow-3m.fits.csv" -80 110 20 0 -p 0.2
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_21_17_41_12-30-0-SRT-Krakow-3m.fits.csv" -60 120 13 0 -l cl -p 0.2
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_21_21_34_12-25-0-SRT-Krakow-3m.fits.csv" -60 120 8 0 -l cl -p 0.2 -r 5
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_21_25_09_12-20-0-SRT-Krakow-3m.fits.csv" -60 120 5 0 -l cl -p 0.2
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_21_28_46_12-15-0-SRT-Krakow-3m.fits.csv" -60 120 8 0 -l cl -p 0.1 -r 2
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_21_32_22_12-10-0-SRT-Krakow-3m.fits.csv" -50 80 10 0 -l cl -p 0.2
python plot_spectrum.py "Krakow_1/pdeverchere-2026_07_12_21_35_48_12-5-0-SRT-Krakow-3m.fits.csv" -50 50 -10 0 -l cl -p 0.2
REM -- A faire 0°