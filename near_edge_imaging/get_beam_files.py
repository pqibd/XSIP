class GetBeamFiles:

    def __init__(self,*,edge_path,flat_path,dark_path,clip):
        flat_files = file_search(flat_path,'*.tif')