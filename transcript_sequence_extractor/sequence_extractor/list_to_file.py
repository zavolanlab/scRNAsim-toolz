def list_to_file(
    to_write_to_file: str,
    filename: str,
) -> None:
    """Creates a file from a list that is input to the function.

    Args:
        to_write_to_file: The list that you want to write to a file.
        filename: The name you want the output file to have (also include the extension of the file while calling the function).

    Returns:
        Nothing, since it outputs a file directly to the working directory
    """
    file = open(filename,'a')   
    for item in to_write_to_file:
        file.write(item)
    file.close()