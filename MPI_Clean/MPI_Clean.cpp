// MPI_Clean.cpp : Definiert den Einstiegspunkt für die Konsolenanwendung.
//

#include "stdafx.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <time.h>
#include <fstream>

struct COO_node {
	int x;
	COO_node* next_node;
};
struct COO_start {
	COO_node* Begin_of_row;
};


int Init(bool *dimension, bool *reading_mode, char* input_dimension, char* input_reading_mode, int argument_count, char* input_COO, bool*coo_mode);
int ErrorOutput(int error, int world_rank);
int SetupGrid(int size, int world_size, int* local_height, int world_rank, bool dimension, bool *unused_processor, int *start_w, int* start_h, int* local_size,int* local_width);
int SetupStartWH(int *start_w,int *start_h,int local_height,int world_rank,int size);
int InitArrays(bool* local_grid, bool* local_buffer, bool*top_halo, bool* bottom_halo, bool* left_halo, bool*right_halo, bool*sending_left, bool*sending_right, int local_size, int local_height, int* true_columns, int size, bool dimension, int * global_true_columns, COO_start* local_grid_coo, bool coo_on);
int fillGridByRandom(bool* grid, int size, double chance);
int fillGridFromData(char* input, int size, int start_w, int start_h, int local_height, bool * local_grid, bool two_dimension, int world_rank, bool COO, COO_start* local_grid_COO);
int SendAndRecieve1D(int world_rank, int size, bool*top_halo, bool* bottom_halo, bool* grid_local, int local_height, int world_size);
int SendAndRecieve2D(int world_rank, int size, bool*top_halo, bool* bottom_halo, bool* left_halo, bool* right_halo, bool* sending_left, bool* sending_right, bool* grid_local, int local_height, int world_size);
int CheckNeighbours1D(int local_size, bool* grid_local, int size, bool* local_buffer, bool*top_halo, bool* bottom_halo, int local_height, int* true_columns, int start_w, int start_h, int* true_diagonal);
int CheckNeighbours2D(int local_size, bool* grid_local, int size, bool* local_buffer, bool*top_halo, bool* bottom_halo, bool*left_halo, bool*right_halo, int local_height, int* true_columns, int start_w, int start_h, int* true_diagonal);
int setBuffer(bool local_value, int count_nachbar, bool*local_buffer, int current_position);
int setColumns(int current_position, bool*local_buffer, int* true_columns, int size, int start_w);
int setDiagonal(int current_position, bool*local_buffer, int start_w, int start_h, int local_height);
int ReviveGrid(int size, bool*grid, int* true_columns, int start_w, int start_h, int local_height, int * true_diagonal);
int ReviveColumn(int size, bool* grid, int column, int start_w, int start_h, int local_height, int * true_diagonal);
int ReviveDiagonale(int size, bool* grid, int start_w, int start_h, int local_height);
int PushBuffer(int size, bool*buffer, bool* local_grid, int world_rank);
int fillGridByRandomCOO(int local_height, COO_start * grid, int sizeOfRow, double chance);
int COO_add_node(int value_x, COO_start* height);



//argv 1: 0=Read in from text(a) or 1= generate Grid(b)
//argv 2: 0=1D or 1=2D
//argv 3: COO, yes or not
//argv 4: Number of Rows or Input File
//argv 5b:  Chance of being alive,
int main(int argc, char** argv)
	//TODO: COO BULLSHIT
{
	// Initialize the MPI environment
	MPI_Init(NULL, NULL);
	// Find out rank, size
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	
	//Declaring all needed variables for init
	bool Read_from_text = false;
	bool two_dimensions = false;
	bool COO = false;
	int error = 0;
	//setting init variables
	error = Init(&two_dimensions, &Read_from_text, argv[2], argv[1], argc, argv[3],&COO);
	//Output Error for init
	if (error > 0) {
		ErrorOutput(error, world_rank);
		return error;
	}
	//Declaring variables needed for SetupGenGrid for both inputs
	int sizeOfRow;
	int local_height=0;
	int local_width = 0;
	bool unused_Processor = false;
	int start_w = 0;
	int start_h = 0;
	int local_size = 0;
	int local_true = 0;
	
	//Checkup for generated random Grid if correct splitting is doable and Splitting
	if (Read_from_text == false) {
		sizeOfRow = atoi(argv[4]);
	}
	else {
		std::string line;
		std::ifstream inputfile(argv[4]);
		if (inputfile.is_open()) {
			std::getline(inputfile,line);
			sizeOfRow = atoi(line.c_str());
		}
		inputfile.close();
	}
	error = SetupGrid(sizeOfRow, world_size, &local_height, world_rank, two_dimensions, &unused_Processor, &start_w, &start_h, &local_size,&local_width);
	//Output error for Checkup generated grid
	if (error > 0) {
		ErrorOutput(error, world_rank);
		return error;
	}
	//variables needed to init arrays
	bool * local_grid=0;
	bool * local_buffer=0;
	bool * top_halo=0;
	bool * bottom_halo=0;
	bool * left_halo=0;
	bool * right_halo=0;
	bool * sending_left=0;
	bool * sending_right=0;
	int * true_columns = 0;
	int * global_true_columns = 0;
	COO_start* starting_array_coo=0;
	//init arrays
	error=InitArrays(local_grid,local_buffer, top_halo, bottom_halo, left_halo, right_halo, sending_left, sending_right,local_size,local_height, true_columns,sizeOfRow,two_dimensions,global_true_columns,starting_array_coo,COO);
	if (error > 0) {
		ErrorOutput(error, world_rank);
		return error;
	}
	//variables needed to fill arrays
	double chance = atof(argv[5]);
	int true_local = 0;
	int global_true_local = 0;
	int true_diagonal = 0;
	int global_true_diagonal = 0;
	//filling the arrays
	if (Read_from_text == false) {		
		//FillGrid
		srand(time(NULL));
		if (COO) {
			error=fillGridByRandomCOO(local_height,starting_array_coo,sizeOfRow,chance);
			if (error > 0) {
				ErrorOutput(error, world_rank);
				return error;
			}
		}
		else {
			fillGridByRandom(local_grid, local_size, chance);
		}
	}
	else {
		error=fillGridFromData(argv[4],sizeOfRow,start_w,start_h,local_height,local_grid,two_dimensions,world_rank,COO,starting_array_coo);
	}
	if (error > 0) {
		ErrorOutput(error, world_rank);
		return error;
	}
	//variables needed to apply rules
	int iterations = 50;
	int CutOff_Grid = (sizeOfRow*sizeOfRow)*0.2;
	int CutOff_Row = sizeOfRow * 0.2;
	//Send Data and apply rules
	for (int i = 0; i < iterations; i++) {
		if (two_dimensions) {
			SendAndRecieve2D(world_rank, sizeOfRow, top_halo, bottom_halo, left_halo, right_halo, sending_left, sending_right, local_grid, local_height, world_size);
			true_local = CheckNeighbours2D(local_size, local_grid, sizeOfRow, local_buffer, top_halo, bottom_halo, left_halo, right_halo, local_height, true_columns, start_w, start_h, &true_diagonal);
		}
		else {
			SendAndRecieve1D(world_rank, sizeOfRow, top_halo, bottom_halo, local_grid, local_height, world_size);
			true_local = CheckNeighbours1D(local_size, local_grid, sizeOfRow, local_buffer, top_halo, bottom_halo, local_height, true_columns, start_w, start_h, &true_diagonal);
		}
		//Apply Reviving rules
		MPI_Allreduce(&true_local, &global_true_local, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		//Revive Grid
		if (global_true_local < CutOff_Grid) {
			ReviveGrid(local_size, local_buffer, true_columns, start_w, start_h, local_height, &true_diagonal);
		}
		//Revive Column
		MPI_Allreduce(true_columns, global_true_columns, sizeOfRow, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		for (int j = 0; j<sizeOfRow; j++) {
			if (global_true_columns[j] < CutOff_Row) {
				ReviveColumn(local_size, local_buffer, j, start_w, start_h, local_height, &true_diagonal);
			}
		}
		//Revivde Diagonale
		MPI_Allreduce(&true_diagonal, &global_true_diagonal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		if (global_true_diagonal < CutOff_Row) {
			ReviveDiagonale(local_size, local_buffer, start_w, start_h, local_height);
		}
		//Push Buffer in array and output alive cells
		PushBuffer(local_size, local_buffer, local_grid, world_rank);

	}
	//end of program
	MPI_Finalize();
    return 0;
}

int Init(bool *dimension, bool *reading_mode, char* input_dimension, char* input_reading_mode, int argument_count, char* input_COO, bool*coo_mode) {
	if (argument_count < 5) {
		return 1;
	}
	if (atoi(input_dimension) == 0) {
		*dimension = false;
	}
	else if(atoi(input_dimension) == 1) {
		*dimension = true;
	}
	else {
		return 2;
	}
	if (atoi(input_reading_mode) == 0) {
		*reading_mode = true;
	}
	else if (atoi(input_reading_mode) == 1) {
		*reading_mode = false;
		if (argument_count < 6) {
			return 1;
		}
	}
	else{
		return 3;
	}
	if (atoi(input_COO) == 0) {
		*coo_mode = false;
	}
	else if(atoi(input_COO)==1){
		*coo_mode = true;
	}
	else {
		return 10;
	}
	return 0;
}
int ErrorOutput(int error, int world_rank) {
	if (world_rank == 0) {
		switch (error) {
			//Init
		case 1: 
			printf("Not enough arguments.\n");
			break;
		case 2: 
			printf("Unknown Dimension mode.\n");
			break;
		case 3: 
			printf("Unknown Input mode.\n");
			break;
		case 10: 
			printf("Unknown COO mode.\n");
			break;
			//Checkup
		case 4:
			printf("Cant split correctly in 2D, not a square amount of PE units\n");
			break;
		case 5:
			printf("Cant split properly Grid is not a multiple of PE units\n");
			break;
		case 6: 
			printf("Failed to allocate memory.\n");
			break;
		case 7:
			printf("File broken, input not a 0 or 1\n");
			break;
		case 8:
			printf("To many Lines in File.\n");
			break;
		case 9:
			printf("Not enough Lines in File.\n");
			break;
		default:
			break;
		}
	}
	return error;
}
int SetupGrid(int size, int world_size, int* local_height, int world_rank, bool dimension, bool *unused_processor, int* start_w,int*start_h, int* local_size, int* local_width) {
		//2D routine
		if (dimension) {
			int root = sqrt(world_size);
			//Check for corret 2Dsplit
			if (world_size%root > 0) {
				return 4;
			}
			if (size%root>0) {
				return 5;
			}
			//seeting correct local height and starting parameters;
			*local_height = size / root;
			*local_width = *local_height;
			SetupStartWH(start_w, start_h, *local_height, world_rank, size);
		}
		//1D routine
		//Too many PE
		if (world_size > size) {
			if (world_rank < size) {
				*local_height = 1;
				if (world_rank == 0) {
					printf("Too many Processing Units, %d Units not used", (world_size-size));
				}
			}
			else {
				*unused_processor = true;
			}
		}
		//good amount of PE
		else{
			*local_height = size / world_size;
			*local_width = size;
			*start_h = *local_height*world_rank;
			if (world_rank==world_size-1){
				*local_height = *local_height + (size%world_size);
			}
		}
		//setting up grid
		if (dimension) {
			*local_size = *local_height* *local_height;
		}
		else {
			*local_size = *local_height*size;
		}
	return 0;
}
int SetupStartWH(int *start_w, int *start_h, int local_height, int world_rank,int size) {
	*start_h = (world_rank / size) * local_height;
	*start_w = (world_rank % size) * local_height;
	return 0;
}
int fillGridByRandom(bool* grid, int size, double chance) {
	int test;
	for (int i = 0; i < size; i++) {
		if (grid[i] == false) {
			test = rand() % 1057;
			if (test < (float)1057 * chance) {
				grid[i] = true;
			}
			else {
				grid[i] = false;
			}
		}
	}
	return 0;
}
int InitArrays(bool* local_grid, bool* local_buffer,bool*top_halo,bool* bottom_halo,bool* left_halo,bool*right_halo,bool*sending_left,bool*sending_right ,int local_size,int local_height, int* true_columns, int size,bool dimension,int * global_true_columns,COO_start* local_grid_coo, bool coo_on) {
	if (coo_on) {
		local_grid_coo = (COO_start*)calloc(local_height, sizeof(COO_start));
	}	
	local_grid = (bool*)calloc(local_size, sizeof(bool));
	local_buffer = (bool*)calloc(local_size, sizeof(bool));
	top_halo = (bool*)calloc(local_height, sizeof(bool));
	bottom_halo = (bool*)calloc(local_height, sizeof(bool));
	
	if (dimension) {
		left_halo = (bool*)calloc(local_height + 2, sizeof(bool));
		right_halo = (bool*)calloc(local_height + 2, sizeof(bool));
		sending_left = (bool*)calloc(local_height + 2, sizeof(bool));
		sending_right = (bool*)calloc(local_height + 2, sizeof(bool));
	}
	true_columns = (int*)calloc(size, sizeof(int));
	global_true_columns = (int*)calloc(size, sizeof(int));
	if (coo_on) {
		if (local_grid_coo == NULL) {
			return 6;
		}
	}	
		if (local_grid == NULL) {
			return 6;
		}
	if ( local_buffer  == NULL) {
		return 6;
	}
	if ( top_halo  == NULL) {
		return 6;
	}
	if ( bottom_halo  == NULL) {
		return 6;
	}
	if (dimension) {
		if (left_halo == NULL) {
			return 6;
		}
		if (right_halo == NULL) {
			return 6;
		}
		if (sending_left == NULL) {
			return 6;
		}
		if (sending_right == NULL) {
			return 6;
		}
	}
	if (true_columns == NULL) {
		return 6;
	}
	if (global_true_columns == NULL) {
		return 6;
	}
	return 0;
}
int fillGridFromData(char* input, int size, int start_w, int start_h,int local_height, bool * local_grid,bool two_dimension, int world_rank,bool COO, COO_start* local_grid_COO) {
	std::ifstream inputfile(input);
	std::string line;
	int read_lines = 0;
	int gridsize = size*size;
	int error=0;
	getline(inputfile, line);
	
		while (getline(inputfile, line)) {
			if (two_dimension) {
			if (((read_lines) % size) >= start_w) {
				if (((read_lines ) % size) < (start_w + local_height)) {
					if (((read_lines ) / size) >= start_h) {
						if (((read_lines ) / size) < (start_h + local_height)) {
							if (line == "1") {
								if (COO) {
									error=COO_add_node((read_lines%size)-start_w,&local_grid_COO[local_height*((read_lines / size) - start_h)]);
									if (error > 0) {
										return 6;
									}
								}
								local_grid[(((read_lines) % size) - start_w) + (local_height*(((read_lines) / size) - start_h))] = true;						
							}
							else if (line == "0") {
								
								local_grid[(((read_lines) % size) - start_w) + (local_height*(((read_lines) / size) - start_h))] = false;
							
							}
							else {
								return 7;
							}
						}
					}
				}
			}
			
		}
		else{
			if (read_lines >= size*(start_h)) {
				if (read_lines < size*(start_h+local_height)) {
					if (line == "1") {
						if (COO) {
							error=COO_add_node((read_lines%size) - start_w, &local_grid_COO[local_height*((read_lines / size) - start_h)]);
							if (error > 0) {
								return 6;
							}
						}
						

							local_grid[read_lines - (size*start_h)] = true;
						
					}
					else if (line == "0") {
						
							local_grid[read_lines - (size*start_h)] = false;
						
					}
					else {
						return 7;
					}
				}
			}
		}
			read_lines = read_lines + 1;
	}
		if (read_lines > gridsize) {
			return 8;
		}
		if (read_lines < gridsize) {
			return 9;
		}
		

	return 0;

}
int SendAndRecieve1D(int world_rank, int size, bool*top_halo, bool* bottom_halo, bool* grid_local, int local_height, int world_size) {
	//send and recieve Halos
	//recv TOP_Halo if not process 0
	if (world_rank != 0) {
		MPI_Recv(top_halo, size, MPI_C_BOOL, world_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	//send top Halo if not bottom process
	if (world_rank != (world_size - 1)) {
		MPI_Send(grid_local + (size*(local_height - 1) * sizeof(bool)), size, MPI_C_BOOL, world_rank + 1, 0, MPI_COMM_WORLD);
		//recieve bottom halo if not last process
		MPI_Recv(bottom_halo, size, MPI_C_BOOL, world_rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	//send bottom halo if not process 0
	if (world_rank != 0) {
		MPI_Send(grid_local, size, MPI_C_BOOL, world_rank - 1, 0, MPI_COMM_WORLD);
	}
	return 0;
}
int SendAndRecieve2D(int world_rank, int size, bool*top_halo, bool* bottom_halo, bool* left_halo, bool* right_halo, bool* sending_left, bool* sending_right, bool* grid_local, int local_height, int world_size) {
	//send and recieve Halos
	int root = sqrt(world_size);
	//recv TOP_Halo if not top process
	if (world_rank >(root - 1)) {
		MPI_Recv(top_halo, local_height, MPI_C_BOOL, world_rank - root, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	//send top Halo and recieve bottom halo if not bottom process
	if (world_rank <world_size - root) {
		MPI_Send(grid_local + (local_height*(local_height - 1)) * sizeof(bool), local_height, MPI_C_BOOL, world_rank + root, 0, MPI_COMM_WORLD);
		//recieve bottom halo if not last process
		MPI_Recv(bottom_halo, local_height, MPI_C_BOOL, world_rank + root, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	//send bottom halo if not top process 
	if (world_rank >(root - 1)) {
		MPI_Send(grid_local, local_height, MPI_C_BOOL, world_rank - root, 1, MPI_COMM_WORLD);
	}
	//calc left and right halos
	for (int i = 1; i <= local_height; i++) {
		sending_left[i] = grid_local[local_height*(i - 1)];
		sending_right[i] = grid_local[(local_height*i) - 1];
	}
	sending_left[0] = top_halo[0];
	sending_right[0] = top_halo[local_height - 1];

	sending_left[local_height + 1] = bottom_halo[0];
	sending_right[local_height + 1] = bottom_halo[local_height - 1];
	//recv right halo
	if ((world_rank % root) <  (root - 1)) {
		MPI_Recv(right_halo, local_height + 2, MPI_C_BOOL, world_rank + 1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	//send right halo
	if ((world_rank % root) > 0) {
		MPI_Send(sending_left, local_height + 2, MPI_C_BOOL, world_rank - 1, 2, MPI_COMM_WORLD);
	}
	//recv left halo
	if ((world_rank % root) > 0) {
		MPI_Recv(left_halo, local_height + 2, MPI_C_BOOL, world_rank - 1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	//send left halo
	if ((world_rank % root) < (root - 1)) {
		MPI_Send(sending_right, local_height + 2, MPI_C_BOOL, world_rank + 1, 3, MPI_COMM_WORLD);
	}
	return 0;
}
int CheckNeighbours1D(int local_size, bool* grid_local, int size, bool* local_buffer, bool*top_halo, bool* bottom_halo, int local_height, int* true_columns, int start_w, int start_h, int* true_diagonal) {
	int local_ones = 0;
	for (int i = 0; i < local_size; i++) {
		int count_nachbar = 0;
		//rechter nachbar
		if (i%size < (size - 1)) {
			if (grid_local[i + 1]) {
				count_nachbar = count_nachbar + 1;
				//printf("R");
			}
			//rechts oben
			if (i < size) {
				if (top_halo[(i%size) + 1]) {
					count_nachbar = count_nachbar + 1;
					//	printf("rO1");
				}
			}
			else {
				if (grid_local[i + 1 - size]) {
					count_nachbar = count_nachbar + 1;
					//	printf("rO2");
				}

			}

			//rechts unten
			if (i >= ((local_height - 1) * size)) {
				if (bottom_halo[(i%size) + 1]) {
					count_nachbar = count_nachbar + 1;
					//	printf("rU1");
				}
			}
			else {
				if (grid_local[i + size + 1]) {
					count_nachbar = count_nachbar + 1;
					//	printf("rU2");
				}
			}
		}

		//linker nachbar
		if (i%size > 0) {
			if (grid_local[i - 1]) {
				count_nachbar = count_nachbar + 1;
				//	printf("L");
			}

			//links oben
			if (i < size) {
				if (top_halo[(i%size) - 1]) {
					count_nachbar = count_nachbar + 1;
					//		printf("lO1");
				}
			}
			else {
				if (grid_local[i - 1 - size]) {
					count_nachbar = count_nachbar + 1;
					//		printf("lO2");
				}

			}

			//links unten
			if (i >= ((local_height - 1) * size)) {
				if (bottom_halo[(i%size) - 1]) {
					count_nachbar = count_nachbar + 1;
					//		printf("lU1");
				}
			}
			else {
				if (grid_local[i + size - 1]) {
					count_nachbar = count_nachbar + 1;
					//	printf("lU2");
				}
			}
		}

		//oberer nachbar
		if (i < size) {
			if (top_halo[i%size]) {
				count_nachbar = count_nachbar + 1;
				//	printf("o1");
			}
		}
		else {
			if (grid_local[i - size]) {
				count_nachbar = count_nachbar + 1;
				//	printf("o2");
			}
		}

		//unterer nachbar
		if (i >= ((local_height - 1) * size)) {
			if (bottom_halo[i%size]) {
				count_nachbar = count_nachbar + 1;
				//	printf("u1");
			}
		}
		else {
			if (grid_local[i + size]) {
				count_nachbar = count_nachbar + 1;
				//	printf("u2");
			}
		}


		//calc buffer
		//	printf("\nI=%d, nachbarn = %d\n",i, count_nachbar);
		local_ones = local_ones + setBuffer(grid_local[i], count_nachbar, local_buffer, i);
		setColumns(i, local_buffer, true_columns, size, start_w);
		*true_diagonal = *true_diagonal + setDiagonal(i, local_buffer, start_w, start_h, local_height);
		
	}
	return local_ones;
}
int CheckNeighbours2D(int local_size, bool* grid_local, int size, bool* local_buffer, bool*top_halo, bool* bottom_halo, bool*left_halo,bool*right_halo ,int local_height, int* true_columns, int start_w, int start_h,int * true_diagonal) {
	int local_ones=0;
	for (int i = 0; i < local_size; i++) {
		int count_nachbar = 0;
		//oberer nachbar
		if (i < local_height) {
			if (top_halo[i]) {
				count_nachbar = count_nachbar + 1;
				//			printf("o1");
			}
		}
		else {
			if (grid_local[i - local_height]) {
				count_nachbar = count_nachbar + 1;
				//			printf("o2");
			}
		}

		//unterer nachbar
		if (i >= local_height*(local_height - 1)) {
			if (bottom_halo[i%local_height]) {
				count_nachbar = count_nachbar + 1;
				//		printf("u1");
			}
		}
		else {
			if (grid_local[i + local_height]) {
				count_nachbar = count_nachbar + 1;
				//			printf("u2");
			}
		}

		//rechter nachbar
		if ((i + 1) % local_height > 0) {
			if (grid_local[i + 1]) {
				count_nachbar = count_nachbar + 1;
				//		printf("r1");
			}
		}
		else {
			if (right_halo[(i / local_height) + 1]) {
				count_nachbar = count_nachbar + 1;
				//		printf("r2");
			}
		}
		//linker nachbar
		if (i%local_height > 0) {
			if (grid_local[i - 1]) {
				count_nachbar = count_nachbar + 1;
				//		printf("l1");
			}
		}
		else
		{
			if (left_halo[(i / local_height) + 1]) {
				count_nachbar = count_nachbar + 1;
				//			printf("l2");
			}
		}
		//oben rechts
		if ((i + 1) % local_height == 0 || i < local_height) {
			if (i < (local_height - 1)) {
				if (top_halo[i + 1]) {
					count_nachbar = count_nachbar + 1;
					//	printf("t1");
				}
			}
			if ((i + 1) % local_height == 0) {
				if (right_halo[i / local_height]) {
					count_nachbar = count_nachbar + 1;
					//			printf("t2");
				}
			}
		}
		else {
			if (grid_local[i - local_height + 1]) {
				count_nachbar = count_nachbar + 1;
				//	printf("t3");
			}
		}
		//oben links
		if (i  % local_height == 0 || i < local_height) {
			if (i%local_height > 0 && i < local_height) {
				if (top_halo[i - 1]) {
					count_nachbar = count_nachbar + 1;
					//		printf("s1");
				}
			}
			if (i  % local_height == 0) {
				if (left_halo[i / local_height]) {
					count_nachbar = count_nachbar + 1;
					//			printf("s2");
				}
			}
		}
		else {
			if (grid_local[i - local_height - 1]) {
				count_nachbar = count_nachbar + 1;
				//	printf("s3");
			}
		}
		// unten rechts
		//			ganz rechts			oder   untere Zeile 
		if ((i + 1) % local_height == 0 || i >= local_height*(local_height - 1)) {
			// wenn  unten 				und			nicht ganz rechts
			if (i >= (local_height*(local_height - 1)) && (i + 1) % local_height > 0) {
				if (bottom_halo[i%local_height + 1]) {
					count_nachbar = count_nachbar + 1;
					//		printf("b1");
				}
			}
			//wenn ganz rechts
			if ((i + 1) % local_height == 0) {
				if (right_halo[(i / local_height) + 2]) {
					count_nachbar = count_nachbar + 1;
					//		printf("b2");
				}
			}
		}
		else {
			if (grid_local[i + local_height + 1]) {
				count_nachbar = count_nachbar + 1;
				//	printf("b3");
			}
		}
		//unten links
		//ganz links				oder		untere Zeile
		if (i  % local_height == 0 || i >= local_height*(local_height - 1)) {
			//untere Zeile								und nicht ganz links
			if (i > (local_height*(local_height - 1)) && i  % local_height > 0) {
				if (bottom_halo[i%local_height - 1]) {
					count_nachbar = count_nachbar + 1;
					//	printf("g1");
				}
			}
			//ganz links
			if (i  % local_height == 0) {
				if (left_halo[(i / local_height) + 2]) {
					count_nachbar = count_nachbar + 1;
					//	printf("g2");
				}
			}
		}
		else {
			if (grid_local[i + local_height - 1]) {
				count_nachbar = count_nachbar + 1;
				//	printf("g3");
			}
		}
		local_ones=local_ones+setBuffer(grid_local[i],count_nachbar,local_buffer,i);
		setColumns(i, local_buffer, true_columns, local_height, start_w);
		*true_diagonal = *true_diagonal + setDiagonal(i,local_buffer,start_w,start_h,local_height);

	}
	return local_ones;
}	
int setBuffer(bool local_value, int count_nachbar, bool*local_buffer, int current_position) {
	if (local_value) {
		if (count_nachbar == 2) {
			local_buffer[current_position] = true;
			return 1;
		}
		else if (count_nachbar == 3) {
			local_buffer[current_position] = true;
			return 1;
		}
		else {
			local_buffer[current_position] = false;
		}
	}
	else {
		if (count_nachbar == 3) {
			local_buffer[current_position] = true;
			return 1;
		}
		else {
			local_buffer[current_position] = false;
		}
	}
	return 0;
}
int setColumns(int current_position, bool*local_buffer, int* true_columns, int size, int start_w) {
	if (local_buffer[current_position]) {
		true_columns[(current_position%size)+start_w] = true_columns[(current_position%size)+start_w] + 1;
	}
	return 0;
}
int setDiagonal(int current_position,bool*local_buffer, int start_w, int start_h, int local_height) {
	if ((current_position%local_height) + start_w == (current_position / local_height) + start_h) {
		if (local_buffer[current_position]) {
			return 1;
		}
	}
	return 0;
}
int ReviveGrid(int size,bool*grid, int* true_columns,int start_w, int start_h, int local_height, int * true_diagonal) {
	int test;
	for (int i = 0; i < size; i++) {
		if (grid[i] == false) {
			test = rand() % 1057;
			if (test < (float)1057 * 0.2) {
				grid[i] = true;
				setColumns(i,grid,true_columns,size,start_w);
				*true_diagonal = *true_diagonal+setDiagonal(i,grid,start_w,start_h,local_height);
			}
			else {
				grid[i] = false;
			}
		}
	}
	return 0;
}
int ReviveColumn(int size, bool* grid, int column, int start_w, int start_h, int local_height, int * true_diagonal) {
	if (column >= start_w&&column<start_w+local_height) {
		int test;
		int init_value = column%local_height;
		for (int i = init_value; i < size; i=i+local_height) {
			if (grid[i] == false) {
				test = rand() % 1057;
				if (test < (float)1057 * 0.2) {
					grid[i] = true;
					*true_diagonal = *true_diagonal + setDiagonal(i, grid, start_w, start_h, local_height);
				}
				else {
					grid[i] = false;
				}
			}
		}
	}
	return 0;
}
int ReviveDiagonale(int size, bool* grid, int start_w, int start_h, int local_height) {
	int test;
	for (int i = 0; i < size; i++) {
		if ((i%local_height) + start_w == (i / local_height) + start_h) {
			if (grid[i] == false) {
				test = rand() % 1057;
				if (test < (float)1057 * 0.2) {
					grid[i] = true;
				}
				else {
					grid[i] = false;
				}
			}
		}
	}
	return 0;
}
int PushBuffer(int size, bool*buffer, bool* local_grid, int world_rank) {
	int end_true = 0;
	for (int i = 0; i < size; i++) {
		local_grid[i] = buffer[i];
		end_true = end_true + 1;
	}
	printf("After Revive %d cells are alive in PU %d.\n", end_true, world_rank);
	return 0;
}
int COO_add_node(int value_x,COO_start* height) {
	COO_node* temp=(COO_node*)calloc(1,sizeof(COO_node));
	temp->x = value_x;
	temp->next_node = height->Begin_of_row;
	height->Begin_of_row = temp;
	if (temp == NULL) {
		return 1;
	}
	return 0;
}