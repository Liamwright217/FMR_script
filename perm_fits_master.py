#Liam Wright 2016
#Small script to plot freq vs permiabilty, find FMR

#To run this script, have this file in the directory with all .dat files and import the code. Then run plot_perm() and you should have a graphical output

import numpy as np
import matplotlib.pyplot as plt
import glob as glob
import os as os




def perm_fitting(Ms,Hk, a_eff, wr, w, gam, Meff):
    
    m = np.multiply(4*np.pi, Ms)
    a = (1+(np.power(a_eff,2)))
    w1 = np.multiply(np.power(wr,2),a) - w
    meff = np.multiply(4*np.pi, Meff)

    
    ur_top_1 = np.multiply((m+Hk),a)
    ur_top_2 = np.multiply(ur_top_1, w1)
    ur_top_3 = ur_top_1 - np.multiply((m-2*Hk),(np.power((np.multiply(a_eff,w)),2)))
    
    ur_top = np.multiply(np.multiply(m,np.power(gam,2)),ur_top_3)

    ur_bot = np.power(w1,2) + np.power(np.multiply(a_eff*gam*w,(m-2*Hk)),2)

    ur = 1 + np.divide(ur_top, ur_bot)
    
    ui_top_1 = np.multiply(a_eff*gam*w,meff)
    ui_top_2 = np.multiply(np.power(np.multiply((meff+Hk),a),2),np.power(gam,2)) + np.power(w,2)

    ui_top = np.multiply(ui_top_1, ui_top_2)
    
    ui_bot = ur_bot = np.power(w1,2) + np.power(np.multiply(a_eff*gam*w,(meff-2*Hk)),2)
    
    ui = np.divide(ui_top, ui_bot)

    return ur, ui




def s21(hard, easy):
    #sets params used in algorithm
    z = 50; c = 170; l = 0.013; t = 100e-9

    #frequency sweep range
    f_easy = easy[:,0];

    #pulls easy real and imag data from .dat file
    s21r_easy = easy[:,5];
    s21i_easy = easy[:,6];

    #creates complex data of real + imag easy
    s21_easy = s21r_easy + np.complex(0,1)*s21i_easy;

    #pulls hard real and imag data from .dat file
    s21r_hard = hard[:,5];
    s21i_hard = hard[:,6];

    #creates complex data of real + imag hard
    s21_hard = s21r_hard + np.complex(0,1)*s21i_hard;
    

    #creates complex conj of easy and hard complex data
    s21_easy_conj = np.conj(s21_easy);
    s21_hard_conj = np.conj(s21_hard);


    #top half of algorithm Kalarickal et al
    top = np.complex(0,1)*np.log(np.divide((s21_easy*s21_hard_conj),(s21_hard*s21_hard_conj)));

    #bottom half of algorithm
    bot = np.log(s21_hard);
    botc = np.conj(bot);

    #permiability
    u = np.divide(top*botc, bot*botc);
    
    return u

    

def s11_s21(hard, easy, z, t, l, c):

    #frequency sweep range
    f_easy = easy[:,0];

    #pulls real and imag data, reflected and transmitted, from .dat file
    s11r_easy = easy[:,1];
    s11i_easy = easy[:,2];
    s21r_easy = easy[:,5];
    s21i_easy = easy[:,6];

    #creates complex data from trans and ref easy data
    s11_easy = s11r_easy + np.complex(0,1)*s11i_easy;
    s21_easy = s21r_easy + np.complex(0,1)*s21i_easy;

    #pulls real and imag data, reflected and transmitted, from .dat file
    s11r_hard = hard[:,1];
    s11i_hard = hard[:,2];
    s21r_hard = hard[:,5];
    s21i_hard = hard[:,6];

    #creats complex data for trans and ref hard data
    s11_hard = s11r_hard + np.complex(0,1)*s11i_hard;
    s21_hard = s21r_hard + np.complex(0,1)*s21i_hard;


    #breaks up Ding et al algorithm
    easy_top = (1 + s11_easy - s21_easy);
    easy_bot = (1 - s11_easy);
    easy_bot_conj = np.conj(easy_bot);
    easy_frac = np.divide(easy_top*easy_bot_conj, easy_bot*easy_bot_conj);
    
    hard_top = (1 + s11_hard - s21_hard);
    hard_bot = (1 - s11_hard);
    hard_bot_conj = np.conj(hard_bot);    
    hard_frac = np.divide(hard_top*hard_bot_conj, hard_bot*hard_bot_conj);
    
    const = z/(np.complex(0,(l*c*t*2*np.pi*4e-7)));

    #permiability
    u = np.multiply(const, (easy_frac - hard_frac));
    
    return u



#function to pull all the filenames into array
def pull_data():

    list_dir = os.listdir('.');
    list_easy = [s for s in list_dir if "easy" in s]

    return list_easy
    
   # for file in glob.glob('*easy*.dat'):
   #     list_easy.append(file)



#plots the perm from data files
def plot_perm():  

    #asks whether the output should be s21 or s11s21 
    a = raw_input('s21 or s11 s21? ')
    list_easy = pull_data()
    
    
    hard_dat = glob.glob('*hard_30mT.dat')

    #iterates through each file and plots the perm vs freq
    for i in range(np.size(list_easy)):
        data_easy = np.genfromtxt(list_easy[i], comments='%', delimiter=',');
        data_hard = np.genfromtxt(hard_dat[0], comments='%', delimiter=',');
        hard = data_hard[0:1400];
        easy = data_easy[0:1400];

        #selects whether s21 or s11s21
        if a == ('s21'):

            u = s21(hard, easy)
        else:
            
            u = s11_s21(hard,easy, z, t, l, c)
    
    
        #splits real and imag
        ureal = np.real(u);
        uimag = np.imag(u);

        #for the legend when plotting, can be commented out but also remove label from plt.plot(...)
        title = list_easy[i].split('_')
        mag = title[3]
        name = title[0]
       
        #plt.plot(easy[:,0], ureal, label=['Real',mag]);
        plt.plot(hard[:,0], uimag, label=['Imag',mag]);
        #plt.show()

    
    plt.title(r'$\mu$ Py81')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel(r'$\mu$')
    plt.legend(loc='upper right')
    plt.grid()
    plt.plot()
    plt.show()


#plot_perm()
#plt.grid()
#plt.show()


