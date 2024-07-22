Converting data into PMO format
================================

The data for the input into the portable microhaoplotype object is often in several files which can either be text file or excel (especially for meta data). pmotools comes with several tools to help convert this input into the various parts of the pmo format and a function to wrap them all together.




Sample Metadata
--------------------------------
Sample metadata can often be in either text files or excel sheets. There two tools that can convert either one into json.

Excel into json

..  code-block:: bash

    pmotools-runner.py excel_meta_to_json_meta -h


Output

    usage: pmotools-runner.py excel_meta_to_json_meta [-h] --file FILE [--sheet SHEET] [--index_col_name INDEX_COL_NAME] --output OUTPUT [--overwrite]

    options:
      -h, --help            show this help message and exit
      --file FILE           Input excel file path
      --sheet SHEET         The sheet to convert, if none provided will default to first sheet
      --index_col_name INDEX_COL_NAME
                            by default output is a list, if an index column name is supplied it will be a dict with this column as index
      --output OUTPUT       Output json file path
      --overwrite           If output file exists, overwrite it


By default it will convert sheet one into json

..  code-block:: bash

    pmotools-runner.py excel_meta_to_json_meta --file metadata.xlsx --output metadata.json --overwrite

To convert a specific sheet supply the sheet name

..  code-block:: bash

    pmotools-runner.py excel_meta_to_json_meta --file metadata.xlsx --output metadata.json --overwrite --sheet Sheet2

By default it will convert the data into a list of objects of the input data, to make it a key dictionary of objects, supply the key/index column

..  code-block:: bash

    pmotools-runner.py excel_meta_to_json_meta --file metadata.xlsx --output metadata.json --overwrite --sheet Sheet2 --index_col_name experiment_id

