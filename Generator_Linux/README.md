# Generator
Info
Perplex based program to produce thermophysical table used in LitMod

** To run this program in Linux based platform you need to install "wine"
. This allows you to run *.exe files.

To install wine run following commands in the terminal

"sudo apt-get install wine-binfmt
"sudo update-binfmts --import /usr/share/binfmts/wine"



Put all the .exe file in your path throught ~/.bashsrc file
 

After you have done as suggested above, you run "Generator_LINUX" from the same directory.
You have to make it executable by changing it permission by running following command.

"chmod 755 Generator_LINUX"

then
"./Generator_LINUX"

After type enter a command prompt will appear and you will be asked certain question which answer and thermophysical property tables will put in the Build_dat_(name of the run you gave). Inside this folder there will be a file called TABLE_LITMOD_use copy this file into the folder your model and assign it a name e.g. 98,97 etc. When you will make models with GUI have to assign these file number to your lithospheric mantle bodies. 

 
