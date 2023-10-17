# Gene Information Retrieval Script

This Python script is designed to retrieve gene information using the [MyGene](https://www.genenames.org/](https://mygene.info)) API and other data sources. It takes a list of gene symbols as input and fetches relevant details for those genes.

## Prerequisites

Before using this script, make sure you have the following software and libraries installed:

- Python 3.x
- Required Python libraries:
  - requests
  - urllib3
  - xml.etree.ElementTree
  - mygene
  - pandas
  - numpy
  - pyodbc
  - mysql.connector

## Usage

1. Clone this repository or download the script to your local machine.

2. Open your terminal or command prompt.

3. Navigate to the directory where the script is located.

4. Execute the script by running the following command: [python script.py input_file.txt output_file_name integer]

  
    - `script.py`: The name of the Python script.
    - `input_file.txt`: A text file containing a list of gene symbols (one symbol per line).
    - `output_libName`: The name of the output file.
    - `sqNumber`: A sequence number used in the output file name (for use in specific pipeline needs).

5. The script will perform gene symbol retrieval and create an output file named `NISD-output_libName-sqNumber.csv` in the same directory.

## Output

The script will generate a CSV file used as input for a pipeline but also a results CSV containing gene information for the provided gene symbols. The columns in the output file may include:

- Query	Not Found
- ID Scores
- Official Symbol
- Official Gene Name
- Alias Input
- Taxonomy ID


## License

This script is open-source and available under the [MIT License](LICENSE).

## Acknowledgments

- This script makes use of the MyGene API to retrieve gene information.
- Thanks to the authors and contributors of the libraries and packages used in this script.


