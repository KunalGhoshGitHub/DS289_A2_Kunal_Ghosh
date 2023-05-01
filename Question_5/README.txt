# PLEASE RUN THE CODE INSIDE AN EMPTY FOLDER TO AVOID ANY CONFLICTS

# Code:
"Q5.cpp"

# Input files:
"Input.txt"
"Value_of_Ns.txt"
"Value_of_rd_s.txt"
"Question_4_Error_Data.csv" (In case this file is deleted there is a copy of this file in the folder "Q_4_Error_Data")
Whenever, an output file in opened and closed it will be displayed in the terminal.
(NOTE: In the input files, lines starting with "//" are comments. Do NOT change the order of the variables in the input files)

# Post-processing file:
"Post_processing.ipynb"

# To compile the code:
make 

# To execute the code:
make execute

# To post process the results:
make make_error_from_Q4_Available:
make post_processing
(NOTE: May require root access)
(NOTE: It will try to open the Jupyter Notebook in your default browser. You have to press "Restart & Run all" in the "Kernel" tab after opening the Jupyter Notebook.)

# To delete the generated files:
make clean
(CAUTION: It will remove all the .csv files from the folder.)

# To delete the generated image files:
make clean_images
(CAUTION: It will remove all the .png files from the folder.)

# To delete the object file:
make clean_object_file

# Make file:
"MakeFile"
(NOTE: "MakeFile" is same for all the questions)

# Header files:
"DS289.h"
(NOTE: "DS289.h" is same for all the questions)

# Image files:
All the image files will be generated in .png format

# Post-processing file:
"Post_processing.ipynb"
(NOTE: Post-processing is done in Jupyter-Notebook using Python 3.9.13)
(NOTE: We tried to explore the case of tolerance 1e-6, all the corresponding results are stored in the folder "TOLERANCE_1e-6")
