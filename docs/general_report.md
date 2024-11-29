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

### Item types

The following types of items are supported:

Type            | Input file type | Description
----------------|-----------------|-------------
`file_table`    | CSV or TSV file | The contents of the file will be rendered as a table.
`image`         | PNG image       | The image will be shown in the report.
`plain_text`    | Plain text file | Presents the contents of the file as-is.
`table`         | JSON file       | Presents a JSON representation of a table as a table in the report.
`file_link`     | Any file        | Presents a hyperlink to the file.

## Ouput files

- `reports/general_html_report/{sample}_{type}.general_report.html`
