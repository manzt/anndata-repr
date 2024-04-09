__all__ = [
    "dataframe_to_table",
]

def dataframe_to_table(dataframe, max_len=10):
    # Initialize the table with the correct class names for styling
    table = """<div class='relative overflow-x-auto my-div'>
    <table class="">
    """

    def make_header(columns):
        # Style the header row according to the provided CSS class names
        header = '<thead class="table-header"><tr>'
        for column in columns:
            header += '<th scope="col" class="column-header">' + column + "</th>"
        header += "</tr></thead>"
        return header

    def make_truncated_data(data, max_len):
        # Function to truncate data if it's longer than 10 characters
        if len(str(data)) > max_len:
            return str(data)[:max_len] + "..."
        return str(data)

    def make_row(row, index, max_len):
        # Style each row according to the provided CSS class names, alternating row color not implemented in CSS
        row_html = f'<tr class="table-row">'
        for value in row:
            row_html += (
                f'<td class="table-cell">{make_truncated_data(value, max_len)}</td>'
            )
        row_html += "</tr>"
        return row_html

    # Construct the table header
    table += make_header(dataframe.columns)

    # Construct each row of the table
    for index, row in enumerate(dataframe.itertuples(index=False), start=1):
        table += make_row(row, index, max_len)

    # Close the table and div tags
    table += "</table></div>"
    return table
