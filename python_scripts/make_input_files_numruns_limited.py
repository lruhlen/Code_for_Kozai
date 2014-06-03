#============================================================
# import packages
#============================================================
import numpy as np


#============================================================
def make_input_files():

    #============================================================
    # Define constants
    #============================================================
    mJup = 1.89813e30
    AU = 1.49597871e13
    secPerYear = 3600.0 * 24.0 * 365.0
    block_one = 'ITMN=     2 ITMX=   300\n JADD=     0 JSUB=     0 NATM=   100  \n Atmx= 1.35E+05 Atmn= 5.00E+03 dTAX= 1.00E-06 L/H = 2.00E-00 \n dLmx= 5.00E-03 dLmn= 1.00E-03 dXmx= 4.00E-02 dXmn= 1.00E-02 \n dPmx= 5.00E-02 dPmn= 1.00E-02 Crad= 1.00E-00 Cwrk= 1.00E-00 \n dZmx= 5.00E-03 dZmn= 1.00E-06 dZdt= 1.00E+01  \n epsP= 1.00E-05 epsR= 1.00E-05 epsL= 1.00E-04 epsT= 1.00E-05 \n SMIN= 1.00E-05       1.00E-05       0.00E+37       1.00E-06 \n SMAX= 0.080000       0.080000        2.00E+37     0.080000  \n dTIM= 1.00E+09 FACT= 1.00E+00 dTMN= 1.00E+04 dTMX= 1.00E+15 \n CHMN= 0.030000 CHMX= 0.120000 XX  = 0.710000 YY  = 0.272930 \n JNOU=      0.2 SIGM= 0.20000 \n'
    block_two = 'ZSTA= 2.00E+33\n'
    block_three = 'VALN= 1.07E-07 QVAL= 3.00E+05 \n       7.10E-01       4.00E-05       3.40E-03       2.72E-05\n       5.02E-04       4.53E-05       5.53E-04       3.84E-04\n       3.32E-06       5.69E-05       9.39E-04       1.26E-03\n       8.25E-03       1.65E-03  \n* \n*  TIME     MFLX   \n   0.00     0.00    \n'

 
    # Declare/set the record numbers and .mod names that correspond to the various planet mass models
    mass_record = {} # approx. masses in units of Mjup, record number = last record of the series (i.e. the oldest, most evolved model at that mass)
    temp = '/Volumes/Data/Work/Research/BodenheimerCode/Code_for_Kozai/inputs/JupMassModels/'
    lower_mass_filename='/Volumes/Data/Work/Research/BodenheimerCode/Code_for_Kozai/outputs/2013/June/June_17_2013/v1/chaindown.mod'

    mass_record[10] = [1,temp+'10MjNF.mod'] # not 100% sure about the record number on this one
    mass_record[3] = [1,temp+'3MjNF.mod'] # not 100% sure about the record number on this one
    mass_record[1] = [1,temp+'1MjNF.mod'] # not 100% sure about the record number on this one
    mass_record[0.5] = [28,lower_mass_filename]
    mass_record[0.4] = [44,lower_mass_filename]
    mass_record[0.3] = [68,lower_mass_filename]
    mass_record[0.2] = [100,lower_mass_filename]
    mass_record[0.1] = [152,lower_mass_filename]
    mass_record[0.075] = [176,lower_mass_filename]
    mass_record[0.055] = [200,lower_mass_filename]
    

    ecc_vals = [0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.75,0.8,0.85,0.9]
    semi_vals = AU * np.array([0.05,0.1,0.15,0.2,0.25,0.3,0.35])
    tKozai_vals = [1e6]
    input_file_number = 1
    input_file_path_base = '/Volumes/Data/Work/Research/BodenheimerCode/Code_for_Kozai/inputs/batch_inputs/'
    #    input_file_path_base = input_file_path_base.split('Code_for_Kozai/')[-1]

    for mass in mass_record:
        start_model_string = mass_record[mass][-1].split('Code_for_Kozai/')[-1]
        start_model_string = ' STARTING MODEL    :  2 {:<}\n'.format(start_model_string)

        moda_string = mass_record[mass][0]
        moda_string = ' MODA={:>6} '.format(moda_string)
        
        for tKozai in tKozai_vals:

            end_time = min(20.0*tKozai, 3e8)
            dtime = tKozai / 20.0
            nmod = int(end_time / dtime)
            nrit = int(end_time / dtime / 40.0)
            #            nrit = 1
            
            tK_string = str('{:.2E}'.format(tKozai))
            tK_string = ' KZPE={:>9} '.format(tK_string)

            nmod_string = str('NMOD={:>6} '.format(nmod))
            nrit_string = str('NRIT={:>6} '.format(nrit))

            for ecc in ecc_vals:
                ecc_string = str('{:.3f}'.format(ecc))
                ecc_string = ' ECCN= {:<8} '.format(ecc_string)

                for semi in semi_vals:
                    semi_string = str('{:.2E}'.format(semi))
                    semi_string = 'SEMI={:>9} '.format(semi_string)
                    
                    input_file_path = input_file_path_base+'run'+str(input_file_number)+'.txt'
                    outfile = open(input_file_path,'w')
                    print input_file_path

                    input_file_path = 'outputs/batch_outputs/run'+str(input_file_number)+'.mod'                    
                    input_file_path_string = ' BINARY OUTPUT FILE:  3 {:<}\n'.format(input_file_path)

                    outfile.write(start_model_string)
                    outfile.write(input_file_path_string)
                    outfile.write(moda_string)
                    outfile.write(nmod_string)
                    outfile.write(nrit_string)
                    outfile.write(block_one)
                    outfile.write(tK_string)
                    outfile.write(block_two)
                    outfile.write(ecc_string)
                    outfile.write(semi_string)
                    outfile.write(block_three)
                    outfile.close()


                    #                  print start_model_string,input_file_path_string,moda_string,nmod_string,nrit_string,block_one,tK_string,block_two,ecc_string,semi_string,block_three 

                    input_file_number = input_file_number + 1
                    
    
    return 0
