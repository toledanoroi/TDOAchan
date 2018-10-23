import scipy.signal.fir_filter_design as spf

# filter1 = spf.firwin(20,[0.01,0.8]);

# need to add
#1. DC cancelation
#2. High Frequency cancellation
#3. find a way to define TOA of a chirp
#4. choose better autucorrelation analog signal


# import pyttsx
# engine = pyttsx.init()
# engine.say('Hello Everyone.')
# engine.say('Welcome to the ultrasonic sound localization system by Roi Toledano and Yarden Avraham.')
# engine.say('Today we will examine our localization ability after we improve the algorithms in this project.')
# engine.say('Hope you will enjoy.')
# engine.runAndWait()



import os

# os.system("mpg321 good.wav")