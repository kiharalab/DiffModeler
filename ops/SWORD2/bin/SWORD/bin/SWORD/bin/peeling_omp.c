#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <getopt.h>
#include <omp.h>
#include <time.h>

#define MAX_ITERATION 64
#define MONOMER '_'

/* Constants for default parameter values */
#define DEFAULT_R2 98
#define DEFAULT_SS2_SIZE 8
#define DEFAULT_MIN_PU_SIZE 30
#define DEFAULT_MAX_PU_SIZE 0
#define DEFAULT_D0 6.0
#define DEFAULT_DELTA 1.5
#define DEFAULT_ONLY_SS 0
#define DEFAULT_PRUNING 0
#define DEFAULT_CUTOFF_PRUNING 0.0
#define DEFAULT_MAX_PU_NUMBER 16

/* Function prototypes */
void cleanup(void);
void parse_pdb(const char *);
void parse_dssp(const char *);
void compute_cumulative_sums(int);
void simple_cutting(void);
void double_cutting(void);
double homogeneity(int, int);
void save_pu(void);
void mutual_information(void);
void measure_coeff(void);
void help(void);
void construct_file_path(const char *, const char *, char *, size_t);
void print_and_free_results(int);

/* Global variables */
double **tab_pcontact = NULL;
double **cum_pcontact = NULL;  // Cumulative sum matrix
double pcontact = 0.0;
int *tab_decoupe = NULL;
int *tab_true_num = NULL;
int *tab_ss2 = NULL;
int iteration = 1;
int nbre_de_coupe = 0;
int max_i1 = 0;
int max_i2 = 0;
int max_j1 = 0;
int max_j2 = 0;
int start = 0;
int end = 0;
int max_start = 0;
int max_end = 0;
int ind = 0;  // Number of residues - 1
double max_coeff_matthews = 0.0;
int size_pu = 0;
int min_size_pu = 0;
int LIMITSIZESS2 = 0;
int MIN_SIZE_PU = 0;
int MAX_SIZE_PU = 0;
double D0 = 0.0;
double DELTA = 0.0;
int ONLYSS2 = 0;
int MAXNUMBEROFPU = 20;
int id_pu = 0;
int nbre_pu = 0;
int new_nb_pu = 0;
int pu[MAX_ITERATION][128][2];
int best_pu = -1;
int current_pu = 0;
double CI = 0.0;
/* Pruning and Homogeneity variables */
int PRUNING = 0;
double CUTOFF_PRUNING = 0.0;
double H_PU1 = 0.0;
double H_PU2 = 0.0;
double H_PU3 = 0.0;
int start_pu = 0;
int end_pu = 0;
/* Global variable for output directory */
char NAME_OUTPUT_DIR[1024] = "";
/* Declare MAXR2 as a global variable */
int MAXR2 = 0;
/* Verbose mode */
int VERBOSE = 0;
/* Global variable for program start time */
double program_start_time = 0.0;

/* Function to print verbose messages */
#define VERBOSE_PRINT(...)      \
    do {                        \
        if (VERBOSE) {          \
            printf(__VA_ARGS__);\
        }                       \
    } while (0)

/* Structure to hold iteration results */
typedef struct {
    double max_cr;
    double min_density;
    double CI;
    double R;
    int num_pu;
    int (*pu_indices)[2]; // array of [start, end] indices for each PU
} IterationResult;

/* Array to store results of each iteration */
IterationResult iteration_results[MAX_ITERATION];

int main(int argc, char *argv[]) {
    /* Record the program start time */
    program_start_time = omp_get_wtime();

    /* Free memory at exit */
    atexit(cleanup);

    /* Variables for command-line arguments */
    char NAME_PDB_FILE[1024] = "";
    char NAME_DSSP_FILE[1024] = "";
    MAXR2 = DEFAULT_R2;
    LIMITSIZESS2 = DEFAULT_SS2_SIZE;
    MIN_SIZE_PU = DEFAULT_MIN_PU_SIZE;
    MAX_SIZE_PU = DEFAULT_MAX_PU_SIZE;
    D0 = DEFAULT_D0;
    DELTA = DEFAULT_DELTA;
    ONLYSS2 = DEFAULT_ONLY_SS;
    PRUNING = DEFAULT_PRUNING;
    CUTOFF_PRUNING = DEFAULT_CUTOFF_PRUNING;
    MAXNUMBEROFPU = DEFAULT_MAX_PU_NUMBER;

    /* ----- Addition Start: Number of CPUs Option ----- */
    /* Initialize NUM_CPU to the maximum number of available threads */
    int NUM_CPU = omp_get_max_threads();
    /* ----- Addition End ----- */

    int opt;
    int option_index = 0;
    static struct option long_options[] = {
        {"pdb-file", required_argument, 0, 'p'},
        {"dssp-file", required_argument, 0, 'd'},
        {"max-r2", required_argument, 0, 'r'},
        {"min-ss-size", required_argument, 0, 's'},
        {"min-pu-size", required_argument, 0, 'l'},
        {"max-pu-size", required_argument, 0, 'm'},
        {"d0-value", required_argument, 0, '0'},
        {"delta-value", required_argument, 0, 't'},
        {"only-ss", required_argument, 0, 'o'},
        {"pruning", required_argument, 0, 'g'},
        {"cutoff-pruning", required_argument, 0, 'c'},
        {"max-pu-number", required_argument, 0, 'n'},
        {"output-directory", required_argument, 0, 'O'},
        {"verbose", no_argument, 0, 'v'},
        {"help", no_argument, 0, 'h'},
        /* ----- Addition Start: Number of CPUs Option ----- */
        {"num-cpu", required_argument, 0, 'C'},
        /* ----- Addition End ----- */
        {0, 0, 0, 0}
    };

    /* Parse command-line arguments */
    while ((opt = getopt_long(argc, argv, "p:d:r:s:l:m:0:t:o:g:c:n:O:vhC:", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'p':
                strncpy(NAME_PDB_FILE, optarg, sizeof(NAME_PDB_FILE) - 1);
                break;
            case 'd':
                strncpy(NAME_DSSP_FILE, optarg, sizeof(NAME_DSSP_FILE) - 1);
                break;
            case 'r':
                MAXR2 = atoi(optarg);
                break;
            case 's':
                LIMITSIZESS2 = atoi(optarg);
                break;
            case 'l':
                MIN_SIZE_PU = atoi(optarg);
                break;
            case 'm':
                MAX_SIZE_PU = atoi(optarg);
                break;
            case '0':
                D0 = atof(optarg);
                break;
            case 't':
                DELTA = atof(optarg);
                break;
            case 'o':
                ONLYSS2 = atoi(optarg);
                break;
            case 'g':
                PRUNING = atoi(optarg);
                break;
            case 'c':
                CUTOFF_PRUNING = atof(optarg);
                break;
            case 'n':
                MAXNUMBEROFPU = atoi(optarg);
                break;
            case 'O':
                strncpy(NAME_OUTPUT_DIR, optarg, sizeof(NAME_OUTPUT_DIR) - 1);
                break;
            case 'v':
                VERBOSE = 1;
                break;
            case 'h':
            default:
                help();
                exit(EXIT_FAILURE);
            case 'C':
                NUM_CPU = atoi(optarg);
                if (NUM_CPU < 1) {
                    fprintf(stderr, "Error: Number of CPUs must be at least 1.\n");
                    exit(EXIT_FAILURE);
                }
                break;
        }
    }

    if (strlen(NAME_PDB_FILE) == 0 || strlen(NAME_DSSP_FILE) == 0) {
        fprintf(stderr, "Error: PDB file and DSSP file are required.\n");
        help();
        exit(EXIT_FAILURE);
    }

    omp_set_num_threads(NUM_CPU);

    VERBOSE_PRINT("Starting the program with the following parameters:\n");
    VERBOSE_PRINT("PDB file: %s\n", NAME_PDB_FILE);
    VERBOSE_PRINT("DSSP file: %s\n", NAME_DSSP_FILE);
    VERBOSE_PRINT("Max R2: %d\n", MAXR2);
    VERBOSE_PRINT("Min SS size: %d\n", LIMITSIZESS2);
    VERBOSE_PRINT("Min PU size: %d\n", MIN_SIZE_PU);
    VERBOSE_PRINT("Max PU size: %d\n", MAX_SIZE_PU);
    VERBOSE_PRINT("D0 value: %.2f\n", D0);
    VERBOSE_PRINT("Delta value: %.2f\n", DELTA);
    VERBOSE_PRINT("Only SS: %d\n", ONLYSS2);
    VERBOSE_PRINT("Pruning: %d\n", PRUNING);
    VERBOSE_PRINT("Cutoff Pruning: %.2f\n", CUTOFF_PRUNING);
    VERBOSE_PRINT("Max PU number: %d\n", MAXNUMBEROFPU);
    VERBOSE_PRINT("Output directory: %s\n", NAME_OUTPUT_DIR);
    VERBOSE_PRINT("Verbose mode: %s\n", VERBOSE ? "Enabled" : "Disabled");
    VERBOSE_PRINT("Number of CPUs (threads): %d\n\n", NUM_CPU);

    /* Parse PDB and DSSP files */
    parse_pdb(NAME_PDB_FILE);
    parse_dssp(NAME_DSSP_FILE);

    /* Compute cumulative sums */
    compute_cumulative_sums(ind);

    nbre_pu = 0; /* At the first iteration, there is one PU */
    iteration = 1;

    if (!VERBOSE) {
        printf("Max_CR Min_Density CI R Num_PUs PU_Delineations\n");
    }
    while (iteration < MAX_ITERATION) {
        /* Record iteration start time */
        double iter_start_time = omp_get_wtime();

        nbre_pu = new_nb_pu;
        new_nb_pu = 0;
        best_pu = -1;

        if (nbre_pu > MAXNUMBEROFPU) {
            VERBOSE_PRINT("Maximum number of PUs reached. Exiting.\n");
            if (VERBOSE) {
                print_and_free_results(iteration - 1);
            }
            exit(EXIT_SUCCESS);
        }

        int continuation = 0;

        for (int x = 0; x <= nbre_pu; x++) {
            current_pu = x;
            if (iteration > 0) {
                start = pu[iteration - 1][x][0];
                end = pu[iteration - 1][x][1];
            }
            size_pu = end - start;
            min_size_pu = MIN_SIZE_PU;
            if (size_pu >= min_size_pu) {
                simple_cutting();
                double_cutting();
            }
        }

        if (best_pu == -1) {
            VERBOSE_PRINT("No further cutting possible. Exiting.\n");
            if (VERBOSE) {
                print_and_free_results(iteration - 1);
            }
            exit(EXIT_SUCCESS);
        }

        save_pu();

        if (MAX_SIZE_PU != 0) {
            continuation = 0;
            for (int x = 0; x <= nbre_pu; x++) {
                current_pu = x;
                if (iteration > 0) {
                    start = pu[iteration - 1][x][0];
                    end = pu[iteration - 1][x][1];
                }
                size_pu = end - start;
                if (size_pu > MAX_SIZE_PU) {
                    continuation = 1;
                }
            }
            if (continuation == 0) {
                VERBOSE_PRINT("All PUs are within the maximum size. Exiting.\n");
                if (VERBOSE) {
                    print_and_free_results(iteration - 1);
                }
                exit(EXIT_SUCCESS);
            }
        }

        if (PRUNING != 0) {
            continuation = 0;
            if (nbre_de_coupe == 1) {
                H_PU1 = homogeneity(max_start, max_i1);
                H_PU2 = homogeneity(max_i2, max_end);
                if (H_PU1 >= CUTOFF_PRUNING || H_PU2 >= CUTOFF_PRUNING) {
                    continuation = 1;
                }
            } else if (nbre_de_coupe == 2) {
                H_PU1 = homogeneity(max_start, max_i1);
                H_PU2 = homogeneity(max_i2, max_j1);
                H_PU3 = homogeneity(max_j2, max_end);
                if (H_PU1 >= CUTOFF_PRUNING || H_PU2 >= CUTOFF_PRUNING || H_PU3 >= CUTOFF_PRUNING) {
                    continuation = 1;
                }
            }
            if (continuation == 0) {
                VERBOSE_PRINT("Pruning criteria not met. Exiting.\n");
                if (VERBOSE) {
                    print_and_free_results(iteration - 1);
                }
                exit(EXIT_SUCCESS);
            }
        }

        measure_coeff();
        mutual_information();

        max_coeff_matthews = 0.0;

        /* Record iteration end time */
        if (VERBOSE) {
            double iter_elapsed = omp_get_wtime() - iter_start_time;
            printf("Iteration %d (Elapsed time: %.1f s)\n", iteration, iter_elapsed);
        }

        iteration++;
    }

    /* Print and free results if verbose mode is enabled */
    if (VERBOSE) {
        print_and_free_results(iteration - 1);
    }

    return 0;
}

void cleanup(void) {
    if (tab_decoupe) free(tab_decoupe);
    if (tab_true_num) free(tab_true_num);
    if (tab_ss2) free(tab_ss2);
    if (tab_pcontact) {
        for (int i = 0; i <= ind; i++) {
            if (tab_pcontact[i]) free(tab_pcontact[i]);
        }
        free(tab_pcontact);
    }
    if (cum_pcontact) {
        for (int i = 0; i <= ind + 1; i++) {
            if (cum_pcontact[i]) free(cum_pcontact[i]);
        }
        free(cum_pcontact);
    }
    if (VERBOSE) {
        /* Calculate and print total runtime */
        double total_time = omp_get_wtime() - program_start_time;
        printf("\nTotal runtime: %.2f seconds\n", total_time);
    }
}

void print_and_free_results(int total_iterations) {
    printf("\nFinal Results:\n");
    printf("--------------------------------------------------------------------------------------\n");
    printf("Iter | Max CR | Min Density |     CI    |     R     | Num PUs | PU Delineations\n");
    printf("--------------------------------------------------------------------------------------\n");
    for (int iter = 1; iter <= total_iterations; iter++) {
        printf("%4d | %-6.2lf | %-11.2lf | %-9.6lf | %-9.6lf | %7d | ",
               iter,
               iteration_results[iter].max_cr,
               iteration_results[iter].min_density,
               iteration_results[iter].CI,
               iteration_results[iter].R,
               iteration_results[iter].num_pu);

        for (int x = 0; x < iteration_results[iter].num_pu; x++) {
            printf("%d-%d", iteration_results[iter].pu_indices[x][0], iteration_results[iter].pu_indices[x][1]);
            if (x < iteration_results[iter].num_pu - 1) {
                printf(", ");
            }
        }
        printf("\n");

        if (iteration_results[iter].pu_indices != NULL) {
            free(iteration_results[iter].pu_indices);
            iteration_results[iter].pu_indices = NULL;
        }
    }
    printf("--------------------------------------------------------------------------------------\n");
}

/* Function to construct file paths */
void construct_file_path(const char *output_dir, const char *filename, char *result_path, size_t result_size) {
    if (strlen(output_dir) > 0) {
        char dir_separator[2];
        if (output_dir[strlen(output_dir) - 1] == '/' || output_dir[strlen(output_dir) - 1] == '\\') {
            dir_separator[0] = '\0';
        } else {
            dir_separator[0] = '/';
            dir_separator[1] = '\0';
        }
        snprintf(result_path, result_size, "%s%s%s", output_dir, dir_separator, filename);
    } else {
        strncpy(result_path, filename, result_size);
        result_path[result_size - 1] = '\0'; // ensure null termination
    }
}

/* Function to parse PDB file and calculate contact matrix */
void parse_pdb(const char *pdb_filename) {
    FILE *pdb_file = fopen(pdb_filename, "r");
    if (pdb_file == NULL) {
        fprintf(stderr, "Error: Cannot open PDB file %s\n", pdb_filename);
        exit(EXIT_FAILURE);
    }

    VERBOSE_PRINT("Parsing PDB file: %s\n", pdb_filename);

    char line[85];
    char chain = 'x';
    char ch = '\0';
    double x = 0, y = 0, z = 0;
    ind = 0;
    int capacity = 1024;  // Initial capacity for dynamic arrays

    /* Dynamic arrays for coordinates */
    double *tab_x = (double *)malloc(capacity * sizeof(double));
    double *tab_y = (double *)malloc(capacity * sizeof(double));
    double *tab_z = (double *)malloc(capacity * sizeof(double));
    char *tab_chain = (char *)malloc(capacity * sizeof(char));

    if (tab_x == NULL || tab_y == NULL || tab_z == NULL || tab_chain == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(EXIT_FAILURE);
    }

    /* Prepare the C-alpha coordinates output file */
    char file_ca_coo[1024];
    construct_file_path(NAME_OUTPUT_DIR, "file_ca_coo.pdb", file_ca_coo, sizeof(file_ca_coo));
    FILE *pfile_ca_coo = fopen(file_ca_coo, "w");
    if (pfile_ca_coo == NULL) {
        fprintf(stderr, "Error: Cannot open output file %s\n", file_ca_coo);
        exit(EXIT_FAILURE);
    }

    while (fgets(line, sizeof(line), pdb_file)) {
        if (strncmp("ATOM", line, 4))
            continue;

        char atname[5] = "";
        sscanf(&line[12], "%4s", atname);
        if (strncmp("CA", atname, 2))
            continue;

        if (chain == 'x') {
            if (line[21] == ' ')
                ch = MONOMER;
            else {
                ch = line[21];
                if (chain == MONOMER)
                    chain = ch;
            }
        }

        sscanf(&line[30], "%lf", &x);
        sscanf(&line[38], "%lf", &y);
        sscanf(&line[46], "%lf", &z);

        if (ind >= capacity) {
            capacity *= 2;
            tab_x = (double *)realloc(tab_x, capacity * sizeof(double));
            tab_y = (double *)realloc(tab_y, capacity * sizeof(double));
            tab_z = (double *)realloc(tab_z, capacity * sizeof(double));
            tab_chain = (char *)realloc(tab_chain, capacity * sizeof(char));
            if (tab_x == NULL || tab_y == NULL || tab_z == NULL || tab_chain == NULL) {
                fprintf(stderr, "Memory allocation failed.\n");
                exit(EXIT_FAILURE);
            }
        }

        tab_chain[ind] = ch;
        tab_x[ind] = x;
        tab_y[ind] = y;
        tab_z[ind] = z;

        /* Write the C-alpha coordinates to the file */
        fprintf(pfile_ca_coo, "%s", line);

        ind++;
    }

    fclose(pdb_file);
    fclose(pfile_ca_coo);
    ind--;  // Adjust index to be zero-based

    /* Allocate and initialize tab_pcontact */
    tab_pcontact = (double **)malloc((ind + 1) * sizeof(double *));
    if (tab_pcontact == NULL) {
        fprintf(stderr, "Memory allocation failed for tab_pcontact.\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i <= ind; i++) {
        tab_pcontact[i] = (double *)calloc((ind + 1), sizeof(double));
        if (tab_pcontact[i] == NULL) {
            fprintf(stderr, "Memory allocation failed for tab_pcontact[%d].\n", i);
            exit(EXIT_FAILURE);
        }
    }

    VERBOSE_PRINT("Calculating contact matrix...\n");

    pcontact = 0.0;

    /* Prepare the probability contact matrix output file */
    char file_proba_contact_mat[1024];
    construct_file_path(NAME_OUTPUT_DIR, "file_proba_contact.mat", file_proba_contact_mat, sizeof(file_proba_contact_mat));
    FILE *pfile_proba_contact = fopen(file_proba_contact_mat, "w");
    if (pfile_proba_contact == NULL) {
        fprintf(stderr, "Error: Cannot open output file %s\n", file_proba_contact_mat);
        exit(EXIT_FAILURE);
    }

    /* Write header for MATLAB-compatible matrix file */
    fprintf(pfile_proba_contact, "# Contact probability matrix\n");

    /* Calculate contact probabilities */
    #pragma omp parallel for schedule(dynamic) reduction(+:pcontact)
    for (int i = 0; i <= ind; i++) {
        for (int j = i; j <= ind; j++) {
            double dx = tab_x[i] - tab_x[j];
            double dy = tab_y[i] - tab_y[j];
            double dz = tab_z[i] - tab_z[j];
            double dt2 = dx * dx + dy * dy + dz * dz;
            double dt = sqrt(dt2);

            double tmp = (dt - D0) / DELTA;
            double p = 1.0 / (1.0 + exp(tmp));

            tab_pcontact[i][j] = p;
            tab_pcontact[j][i] = p; // Since the matrix is symmetric

            if (i != j) {
                pcontact += 2 * p;
            } else {
                pcontact += p;
            }
        }
    }
    /* Write the probability contact matrix in matrix format */
    for (int i = 0; i <= ind; i++) {
        for (int j = 0; j <= ind; j++) {
            fprintf(pfile_proba_contact, "%7.5f ", tab_pcontact[i][j]);
        }
        fprintf(pfile_proba_contact, "\n");  // New line at the end of each row
    }

    fclose(pfile_proba_contact);

    /* Initialize variables for PUs */
    start = 0;
    end = ind;
    pu[0][0][0] = start;
    pu[0][0][1] = end;
    iteration = 1;
    best_pu = -1;
    new_nb_pu = 0;
    nbre_pu = 0;

    /* Free coordinate arrays */
    if (tab_x) free(tab_x);
    if (tab_y) free(tab_y);
    if (tab_z) free(tab_z);
    if (tab_chain) free(tab_chain);
}

/* Function to compute cumulative sums */
void compute_cumulative_sums(int size) {
    VERBOSE_PRINT("Computing cumulative sums...\n");

    /* Allocate memory for cum_pcontact */
    //VERBOSE_PRINT("Coucou %ld\n", (size + 2) * sizeof(double *));
    cum_pcontact = (double **)malloc((size + 2) * sizeof(double *));
    if (cum_pcontact == NULL) {
        fprintf(stderr, "Memory allocation failed for cum_pcontact.\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i <= size + 1; i++) {
        cum_pcontact[i] = (double *)calloc(size + 2, sizeof(double));
        if (cum_pcontact[i] == NULL) {
            fprintf(stderr, "Memory allocation failed for cum_pcontact[%d].\n", i);
            exit(EXIT_FAILURE);
        }
    }

    /* Compute the cumulative sums */
    for (int i = 1; i <= size + 1; i++) {
        for (int j = 1; j <= size + 1; j++) {
            cum_pcontact[i][j] = tab_pcontact[i - 1][j - 1] +
                                 cum_pcontact[i - 1][j] +
                                 cum_pcontact[i][j - 1] -
                                 cum_pcontact[i - 1][j - 1];
        }
    }
}

/* Function to get sum over a rectangle using cumulative sums */
double get_rectangle_sum(int row_start, int col_start, int row_end, int col_end) {
    // Adjust indices for cumulative sum matrix
    row_start++;
    col_start++;
    row_end++;
    col_end++;

    double sum = cum_pcontact[row_end][col_end] -
                 cum_pcontact[row_start - 1][col_end] -
                 cum_pcontact[row_end][col_start - 1] +
                 cum_pcontact[row_start - 1][col_start - 1];
    return sum;
}

/* Function to parse DSSP file */
void parse_dssp(const char *dssp_filename) {
    FILE *dssp_file = fopen(dssp_filename, "r");
    if (dssp_file == NULL) {
        fprintf(stderr, "Error: Cannot open DSSP file %s\n", dssp_filename);
        exit(EXIT_FAILURE);
    }

    VERBOSE_PRINT("Parsing DSSP file: %s\n", dssp_filename);

    char line[200];
    int index2 = 0;
    int capacity = ind + 1;

    /* Allocate memory for tab_ss2 and tab_true_num */
    tab_ss2 = (int *)calloc(capacity, sizeof(int));
    tab_true_num = (int *)calloc(capacity, sizeof(int));
    if (tab_ss2 == NULL || tab_true_num == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(EXIT_FAILURE);
    }

    while (fgets(line, sizeof(line), dssp_file)) {
        if (isspace(line[126]) || line[126] == '-')
            continue;

        /* Adjust the line as in the old code */
        line[6] = ' ';
        line[11] = ' ';

        int numaa2;
        sscanf(&line[6], "%d", &numaa2);

        /* Read and check the amino acid type */
        char aatype[2];
        aatype[0] = line[13];
        aatype[1] = '\0';

        if (aatype[0] == '!')
            continue;

        if (index2 >= capacity) {
            capacity *= 2;
            tab_ss2 = (int *)realloc(tab_ss2, capacity * sizeof(int));
            tab_true_num = (int *)realloc(tab_true_num, capacity * sizeof(int));
            if (tab_ss2 == NULL || tab_true_num == NULL) {
                fprintf(stderr, "Memory allocation failed.\n");
                exit(EXIT_FAILURE);
            }
        }

        tab_true_num[index2] = numaa2;

        char ss2 = line[16];

        if (isspace(ss2))
            tab_ss2[index2] = 0;
        else if (ss2 == 'H' || ss2 == 'G')
            tab_ss2[index2] = 1;
        else if (ss2 == 'E' || ss2 == 'B')
            tab_ss2[index2] = 2;
        else
            tab_ss2[index2] = 0;

        index2++;
    }

    index2--;
    fclose(dssp_file);

    /* Initialize tab_decoupe with the larger of ind and index2 */
    int max_residues = (ind > index2) ? ind : index2;
    tab_decoupe = (int *)malloc((max_residues + 1) * sizeof(int));
    if (tab_decoupe == NULL) {
        fprintf(stderr, "Memory allocation failed for tab_decoupe.\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i <= max_residues; i++)
        tab_decoupe[i] = 1;

    /* Process secondary structure segments */
    int start1 = 0;
    int sizess2;
    int startss2 = 0, endss2;
    int tab_u_ss2[1024];
    int j = 0;

    for (int i = 0; i <= index2; i++) {
        if (tab_ss2[i] != 0 && start1 == 0) {
            start1 = tab_ss2[i];
            startss2 = i;
        }
        if (tab_ss2[i] != start1 && tab_ss2[i] != 0) {
            start1 = tab_ss2[i];
            endss2 = i - 1;
            sizess2 = endss2 - startss2;
            if (sizess2 <= LIMITSIZESS2 && sizess2 != 0) {
                tab_u_ss2[j++] = startss2;
                tab_u_ss2[j++] = endss2;
            }
        }
        if (tab_ss2[i] == 0 && start1 != 0) {
            start1 = 0;
            endss2 = i - 1;
            sizess2 = endss2 - startss2;
            if (sizess2 <= LIMITSIZESS2 && sizess2 != 0) {
                tab_u_ss2[j++] = startss2;
                tab_u_ss2[j++] = endss2;
            }
        }
    }

    for (int i = 0; i < j; i += 2)
        for (int k = tab_u_ss2[i]; k < tab_u_ss2[i + 1]; k++)
            tab_decoupe[k] = 0;

    if (ONLYSS2 == 1) {
        for (int i = 0; i <= ind; i++)
            for (int j = 0; j <= ind; j++)
                if (tab_ss2[i] == 0 || tab_ss2[j] == 0)
                    tab_pcontact[i][j] /= 10.0;
    }
}

/* Simple cutting function */
void simple_cutting(void) {
    int min_seg_size = MIN_SIZE_PU;
    for (int i = start; i < end; i++) {
        int i_p = i + 1;
        if (tab_decoupe[i_p] == 0)
            continue;
        int size_pu1 = i - start;
        int size_pu2 = end - i;
        if (size_pu1 >= min_seg_size && size_pu2 >= min_seg_size) {
            int start_ju = start;
            int i1 = i;
            int i2 = i + 1;
            int end_ju = end;
            double a = get_rectangle_sum(start_ju, start_ju, i1, i1);
            double b = get_rectangle_sum(i2, i2, end_ju, end_ju);
            double c = get_rectangle_sum(start_ju, i2, i1, end_ju);

            double denom = (a + c) * (b + c);
            if (denom == 0.0) continue; // Avoid division by zero

            double coeff_matthews = (a * b - c * c) / denom;
            if (coeff_matthews > max_coeff_matthews) {
                max_start = start_ju;
                max_i1 = i1;
                max_i2 = i2;
                max_end = end_ju;
                max_coeff_matthews = coeff_matthews;
                nbre_de_coupe = 1;
                best_pu = current_pu;
            }
        }
    }
}

/* Modified Double cutting function */
void double_cutting(void) {
    int min_seg_size = MIN_SIZE_PU;
    int coo_max_i = end - min_seg_size;
    int coo_max_j = end - (int)(min_seg_size / 2);
    int coo_min_i = start + min_seg_size - 1;

    double local_max_coeff_matthews = max_coeff_matthews;

    int completed_iterations = 0;

    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic)
        for (int i = coo_min_i; i < coo_max_i; i++) {
            if (tab_decoupe[i + 1] == 0)
                continue;
            int coo_min_j = i + min_seg_size;
            if (coo_min_j + min_seg_size < coo_max_j) {
                for (int j = coo_min_j; j <= coo_max_j; j++) {
                    if (tab_decoupe[j + 1] == 0)
                        continue;
                    int size_pu2 = j - i;
                    if (size_pu2 >= min_seg_size) {
                        int start_ju = start;
                        int i1 = i;
                        int i2 = i + 1;
                        int j1 = j;
                        int j2 = j + 1;
                        int end_ju = end;

                        // Compute 'a'
                        double a = get_rectangle_sum(i2, i2, j1, j1);

                        // Compute 'b'
                        double b1 = get_rectangle_sum(start_ju, start_ju, i1, i1);
                        double b2 = get_rectangle_sum(j2, j2, end_ju, end_ju);
                        double b3 = get_rectangle_sum(start_ju, j2, i1, end_ju);
                        double b = b1 + b2 + 2 * b3;

                        // Compute 'c'
                        double c1 = get_rectangle_sum(i2, start_ju, j1, i1);
                        double c2 = get_rectangle_sum(j2, i2, end_ju, j1);
                        double c = c1 + c2;

                        double denom = (a + c) * (b + c);
                        if (denom == 0.0) continue; // Avoid division by zero

                        double coeff_matthews = (a * b - c * c) / denom;

                        #pragma omp critical
                        {
                            if (coeff_matthews > local_max_coeff_matthews) {
                                local_max_coeff_matthews = coeff_matthews;
                                max_start = start_ju;
                                max_i1 = i1;
                                max_i2 = i2;
                                max_j1 = j1;
                                max_j2 = j2;
                                max_end = end_ju;
                                max_coeff_matthews = coeff_matthews;
                                nbre_de_coupe = 2;
                                best_pu = current_pu;
                            }
                        }
                    }
                }
            }

            #pragma omp atomic
            completed_iterations++;
        }
        // Wait for all threads to finish
        #pragma omp barrier
    }
}

/* Function to save PUs */
void save_pu(void) {
    if (nbre_de_coupe == 1) {
        pu[iteration][new_nb_pu][0] = max_start;
        pu[iteration][new_nb_pu][1] = max_i1;
        new_nb_pu++;
        pu[iteration][new_nb_pu][0] = max_i2;
        pu[iteration][new_nb_pu][1] = max_end;
    } else if (nbre_de_coupe == 2) {
        pu[iteration][new_nb_pu][0] = max_start;
        pu[iteration][new_nb_pu][1] = max_i1;
        new_nb_pu++;
        pu[iteration][new_nb_pu][0] = max_i2;
        pu[iteration][new_nb_pu][1] = max_j1;
        new_nb_pu++;
        pu[iteration][new_nb_pu][0] = max_j2;
        pu[iteration][new_nb_pu][1] = max_end;
    }

    if (iteration > 1) {
        for (int x = 0; x <= nbre_pu; x++) {
            if (x != best_pu) {
                new_nb_pu++;
                pu[iteration][new_nb_pu][0] = pu[iteration - 1][x][0];
                pu[iteration][new_nb_pu][1] = pu[iteration - 1][x][1];
            }
        }
    }
}

/* Function to calculate homogeneity */
double homogeneity(int start_pu, int end_pu) {
    double pcontact1 = 0.0;
    double pcontact2 = 0.0;
    double threshold = 0.5;
    for (int k = start_pu; k <= end_pu; k++)
        for (int l = start_pu; l <= end_pu; l++)
            if (tab_pcontact[k][l] > threshold) {
                pcontact1 += tab_pcontact[k][l];
                if (abs(k - l) < 6)
                    pcontact2 += tab_pcontact[k][l];
            }

    double H_1 = 0.0, H_2 = 0.0;
    for (int i = start_pu; i < end_pu; i++) {
        for (int j = start_pu; j < end_pu; j++) {
            if (tab_pcontact[i][j] < 0.0001)
                continue;
            if (tab_pcontact[i][j] > threshold) {
                double pnormalize1 = tab_pcontact[i][j] / pcontact1;
                H_1 += pnormalize1 * log(pnormalize1);
                if (abs(i - j) < 6) {
                    double pnormalize2 = tab_pcontact[i][j] / pcontact2;
                    H_2 += pnormalize2 * log(pnormalize2);
                }
            }
        }
    }

    double Neq1 = exp(-H_1);
    double Neq2 = exp(-H_2);
    int n = end_pu - start_pu;
    double I = (Neq1 - Neq2) / n;
    return I;
}

/* Function for mutual information */
void mutual_information(void) {
    double prob_zone[64][64] = {{0.0}};
    double sprob_zone[64] = {0.0};
    double sprob_tot = 0.0;
    char file_matrix_pu_contact[1024];
    char file_pu_delineation[1024];

    construct_file_path(NAME_OUTPUT_DIR, "file_matrix_pu_contact.mtx", file_matrix_pu_contact, sizeof(file_matrix_pu_contact));
    construct_file_path(NAME_OUTPUT_DIR, "file_pu_delineation.mtx", file_pu_delineation, sizeof(file_pu_delineation));

    FILE *pfile_pu_contact_ie = fopen(file_matrix_pu_contact, "w");
    if (pfile_pu_contact_ie == NULL) {
        fprintf(stderr, "Error: Cannot open output file %s\n", file_matrix_pu_contact);
        exit(EXIT_FAILURE);
    }
    FILE *pfile_pu_delineation = fopen(file_pu_delineation, "w");
    if (pfile_pu_delineation == NULL) {
        fprintf(stderr, "Error: Cannot open output file %s\n", file_pu_delineation);
        exit(EXIT_FAILURE);
    }

    for (int x = 0; x <= new_nb_pu; x++) {
        int x1 = pu[iteration][x][0];
        int x2 = pu[iteration][x][1];
        for (int y = 0; y <= new_nb_pu; y++) {
            int y1 = pu[iteration][y][0];
            int y2 = pu[iteration][y][1];
            prob_zone[x][y] = 0.0;
            for (int k1 = x1; k1 <= x2; k1++)
                for (int k2 = y1; k2 <= y2; k2++)
                    prob_zone[x][y] += tab_pcontact[k1][k2];
            fprintf(pfile_pu_contact_ie, "%d %d %f\n", x, y, prob_zone[x][y]);
        }
        fprintf(pfile_pu_delineation, "%d %d %d\n", x, x1, x2);
    }
    fclose(pfile_pu_contact_ie);
    fclose(pfile_pu_delineation);

    for (int x = 0; x <= new_nb_pu; x++) {
        sprob_zone[x] = 0.0;
        for (int y = 0; y <= new_nb_pu; y++)
            sprob_zone[x] += prob_zone[x][y];
        sprob_tot += sprob_zone[x];
    }

    for (int x = 0; x <= new_nb_pu; x++) {
        for (int y = 0; y <= new_nb_pu; y++)
            prob_zone[x][y] /= sprob_tot;
        sprob_zone[x] /= sprob_tot;
    }

    double entropy = 0.0;
    for (int x = 0; x <= new_nb_pu; x++)
        for (int y = 0; y <= new_nb_pu; y++)
            if (prob_zone[x][y] > 0.00001 && sprob_zone[x] > 0.00001 && sprob_zone[y] > 0.00001)
                entropy += prob_zone[x][y] * log(prob_zone[x][y] / (sprob_zone[x] * sprob_zone[y]));

    CI = 100 * sqrt(1 - exp(-2 * entropy));
    double R = 100 * (1 - exp(-2 * entropy));

    iteration_results[iteration].CI = CI;
    iteration_results[iteration].R = R;
    iteration_results[iteration].num_pu = new_nb_pu + 1;

    // Allocate memory for pu_indices
    iteration_results[iteration].pu_indices = malloc((new_nb_pu + 1) * sizeof(int[2]));
    if (iteration_results[iteration].pu_indices == NULL) {
        fprintf(stderr, "Memory allocation failed for pu_indices.\n");
        exit(EXIT_FAILURE);
    }

    for (int x = 0; x <= new_nb_pu; x++) {
        iteration_results[iteration].pu_indices[x][0] = tab_true_num[pu[iteration][x][0]];
        iteration_results[iteration].pu_indices[x][1] = tab_true_num[pu[iteration][x][1]];
    }

    if (!VERBOSE) {
        printf("%lf %lf ", CI, R);
        printf("%d ", new_nb_pu + 1);
        for (int x = 0; x <= new_nb_pu; x++)
            printf("%d %d ", tab_true_num[pu[iteration][x][0]], tab_true_num[pu[iteration][x][1]]);
        printf("\n");
    }

    if (CI > MAXR2) {
        VERBOSE_PRINT("\nCI (%lf) exceeds MAXR2 (%d). Exiting.\n", CI, MAXR2);
        if (VERBOSE) {
            print_and_free_results(iteration);
        }
        exit(EXIT_SUCCESS);
    }
}


/* Function to measure coefficients */
void measure_coeff(void) {
    double prob_zone[64][64] = {{0.0}};
    double min_density = 1e6;
    double max_cr = -1e6;

    for (int x = 0; x <= new_nb_pu; x++) {
        int x1 = pu[iteration][x][0];
        int x2 = pu[iteration][x][1];
        prob_zone[x][x] = 0.0;

        int size_PU1 = x2 - x1 + 1;
        for (int k1 = x1; k1 <= x2; k1++)
            for (int k2 = x1; k2 <= x2; k2++)
                prob_zone[x][x] += tab_pcontact[k1][k2];

        double prob_zone_internal_PU1 = prob_zone[x][x];
        double density_PU1 = prob_zone_internal_PU1 / size_PU1;

        if (min_density > density_PU1)
            min_density = density_PU1;

        for (int y = 0; y <= new_nb_pu; y++) {
            if (x == y)
                continue;
            int y1 = pu[iteration][y][0];
            int y2 = pu[iteration][y][1];
            prob_zone[x][y] = 0.0;

            int size_PU2 = y2 - y1 + 1;
            for (int k1 = y1; k1 <= y2; k1++)
                for (int k2 = y1; k2 <= y2; k2++)
                    prob_zone[y][y] += tab_pcontact[k1][k2];

            double prob_zone_internal_PU2 = prob_zone[y][y];

            for (int k1 = x1; k1 <= x2; k1++)
                for (int k2 = y1; k2 <= y2; k2++)
                    prob_zone[x][y] += tab_pcontact[k1][k2];

            double prob_zone_external = prob_zone[x][y];

            double pdp_criterion = (prob_zone_internal_PU1 + prob_zone_internal_PU2 + prob_zone_external) / (size_PU1 + size_PU2);
            double nnc = prob_zone_external / (pow(size_PU1, 0.43) * pow(size_PU2, 0.43));
            double cr = nnc / pdp_criterion;

            if (max_cr < cr)
                max_cr = cr;
        }
    }

    iteration_results[iteration].max_cr = max_cr;
    iteration_results[iteration].min_density = min_density;

    if (!VERBOSE) {
        printf("%-5.2lf %-5.2lf ", max_cr, min_density);
    }
}

/* Help function */
void help(void) {
    printf("Peeling: a tool for cutting a protein into compact subfragments\n");
    printf("by Jean-Christophe Gelly\n\n");
    printf("You agree that you will not distribute this software or any derivative works of this software to any party whatsoever, modified or unmodified.\n\n");
    printf("Publications:\n");
    printf("\t2006; Gelly, J. C.; de Brevern, A. G.; Hazout, S.\n");
    printf("\tProtein Peeling: an approach for splitting a 3D protein structure into compact fragments\n");
    printf("\tBIOINFORMATICS; 22; (2) 129-33; [pmid:16301202]\n\n");
    printf("\t2006; Gelly, J. C.; Etchebest, C.; Hazout, S.; de Brevern, A. G.\n");
    printf("\tProtein Peeling 2: a web server to convert protein structures into series of protein units\n");
    printf("\tNUCLEIC ACIDS RES; 34; (Web Server issue) W75-8; [pmid:16845113]\n\n");
    printf("\t2011; Gelly, J. C.; de Brevern, A. G.\n");
    printf("\tProtein Peeling 3D: new tools for analyzing protein structures\n");
    printf("\tBIOINFORMATICS; 27; (2) 132-133\n\n");
    printf("Usage:\n");
    printf("\tpeeling [options]\n\n");
    printf("Options:\n");
    printf("\t--pdb-file, -p <file>\t\tSet input PDB file (required)\n");
    printf("\t--dssp-file, -d <file>\t\tSet input DSSP file (required)\n");
    printf("\t--max-r2, -r <value>\t\tSet maximum R2 value (default: %d)\n", DEFAULT_R2);
    printf("\t--min-ss-size, -s <value>\tSet minimum size of secondary structure segment that can be cut (default: %d)\n", DEFAULT_SS2_SIZE);
    printf("\t--min-pu-size, -l <value>\tSet minimum size of Protein Unit (PU) (default: %d)\n", DEFAULT_MIN_PU_SIZE);
    printf("\t--max-pu-size, -m <value>\tSet maximum size of PU (default: %d)\n", DEFAULT_MAX_PU_SIZE);
    printf("\t--d0-value, -0 <value>\t\tSet D0 value in logistic function (default: %.1f)\n", DEFAULT_D0);
    printf("\t--delta-value, -t <value>\tSet Delta value in logistic function (default: %.1f)\n", DEFAULT_DELTA);
    printf("\t--only-ss, -o <0|1>\t\tCompute only C-alpha included in regular secondary structure (default: %d)\n", DEFAULT_ONLY_SS);
    printf("\t--pruning, -g <0|1>\t\tEnable pruning (default: %d)\n", DEFAULT_PRUNING);
    printf("\t--cutoff-pruning, -c <value>\tSet cutoff pruning value (default: %.1f)\n", DEFAULT_CUTOFF_PRUNING);
    printf("\t--max-pu-number, -n <value>\tSet maximum number of PUs (default: %d)\n", DEFAULT_MAX_PU_NUMBER);
    printf("\t--output-directory, -O <dir>\tSet output directory\n");
    printf("\t--num-cpu, -C <value>\t\tSet number of CPUs (threads) to use (default: maximum available)\n");
    printf("\t--verbose, -v\t\t\tEnable verbose mode\n");
    printf("\t--help, -h\t\t\tDisplay this help and exit\n\n");
    printf("Example:\n");
    printf("\tpeeling -p input.pdb -d input.dssp -r 98 -s 8 -l 30 -m 0 -0 6.0 -t 1.5 -o 0 -g 0 -c 0 -n 16 -C 4\n\n"); // Updated example with -C
    printf("DISCLAIMER:\n");
    printf("THIS SOFTWARE IS PROVIDED \"AS IS\" AND ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.\n");
    printf("IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE.\n");
}