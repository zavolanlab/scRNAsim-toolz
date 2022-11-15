def list_to_file(to_write_to_file,filename):
    file = open(filename,'a')   
    for item in to_write_to_file:
        file.write(item)
    file.close()