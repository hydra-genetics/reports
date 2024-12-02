# General report

## Configuration

The contents of the general report is defined by the configuration.
This is defined under the configuration key `general_report` and should be a path to a yaml file that lists the items that should be presented in the report.
The yaml file should define an array of object that can have the following properties:

Property      | Required | Description
--------------|----------|-------------
`name`        | yes      | The name of the item. Will the header for this item in the report.
`description` | yes      | A description of the results. Will follow immediately after the header.
`type`        | yes      | The result type. Read more below.
`input`       | yes      | A relative path to the file defining the results.
`nav_header`  | yes      | The header in the navigation bar that the item should be a member of.
`sections`    | no       | Only used for type `multiqc`. Sections in the MultiQC report that should be extracted. Currently only `table` is supported.

The `input` field can contain `sample` and `type` wildcards.

### Item types

The following types of items are supported:

Type            | Input file type   | Description
----------------|-------------------|-------------
`file_table`    | CSV or TSV file   | The contents of the file will be rendered as a table.
`image`         | PNG image         | The image will be shown in the report.
`plain_text`    | Plain text file   | Presents the contents of the file as-is.
`table`         | JSON file         | Presents a JSON representation of a table as a table in the report. See format description below.
`file_link`     | Any file          | Presents a hyperlink to the file.
`multiqc`       | MultiQC JSON file | Extracts and presents results for this sample from a MultiQC report.

#### JSON tables

The JSON for the tables should be represented as an array of objects where each object represents a row.
For example, this table

column1 | column2
--------|---------
row1    | value1
row2    | value2
row3    | value3

would be represented as

```json
[
    {
        "column1": "row1",
        "column2": "value1",
    },
    {
        "column1": "row2",
        "column2": "value2",
    },
    {
        "column1": "row3",
        "column2": "value3",
    }
]
```

## Example configuration

```yaml
files:
  - name: file_table
    type: file_table
    description: Comma or tab separated file
    input: path/to/file_table.csv
    nav_header: File types
  - name: plain_text
    type: plain_text
    description: plain text file
    input: path/to/plain_text.txt
    nav_header: File types
  - name: file_table (txt Tab-separated)
    type: file_table
    description: A tab separated txt file
    input: path/to/table.txt
    nav_header: File types
  - name: table
    type: table
    input: path/to/table.json
    description: Json file with json representation of a table
    nav_header: File types
  - name: file_link
    type: file_link
    input: reports/cnv_html_report/{sample}_{type}.pathology.cnv_report.html
    description: File link to other results
    nav_header: File types
  - name: image
    type: image
    input: path/to/image.png
    description: Image file
    nav_header: File types
  - name: MultiQC
    type: multiqc
    input: qc/multiqc/multiqc_DNA_data/multiqc_data.json
    sections: ["table"]
    description: multiqc general stats for this sample
    nav_header: File types
```

## Ouput files

- `reports/general_html_report/{sample}_{type}.general_report.html`
