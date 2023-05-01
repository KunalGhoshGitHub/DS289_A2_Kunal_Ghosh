.PHONY: start generate clean clean_object_file execute

all: start generate

object_file = Q5.out
code = Q5.cpp
post_processing_file = Post_processing.ipynb
start: 
	@echo "We are running the make file.... \n"

generate: ${code}
	@echo "Compiling the code.... \n"
	g++ -std=c++11 ${code} -o ${object_file}

execute: ${object_file}
	@echo "Executing the code: \n"
	./${object_file}
	
post_processing: ${post_processing_file}
	@echo "We are opening the Jupyter Notebook in your default web browser.... \n"
	jupyter-notebook Post_processing.ipynb

clean:
	@echo "Removing the generated text files.... \n"
	rm *.csv
	
clean_images:
	@echo "Removing the generated image files.... \n"
	rm *.png
	
clean_object_file:
	@echo "Removing the object file.... \n"
	rm ${object_file}
	
make_error_from_Q4_Available:
	@echo "Making Error data from Question 4 available.... \n"
	cp Q_4_Error_Data/Question_4_Error_Data.csv Question_4_Error_Data.csv
