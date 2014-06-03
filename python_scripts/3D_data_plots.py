def plot_mass_period_radius(data): # The input (data) is a pandas dataframe

    files = data.Path.values.tolist() # File paths to each simulation result
    masses = data.Zmass.values.tolist() # Planet masses
    semi = data.SEMI.values.tolist() # Semi-major axis of planets' orbits
    star_masses = data.ZSTA.values.tolist()  # Mass of central stars
    
    
    return 0

def check_num_sims_in_dir(filepath):
    if len(os.listdir(filepath)) <= 2:
        #        print filepath+' has no (or few) converged models, so it is not getting plotted.'
        return False
    if not os.path.exists(filepath+"/atmos_data.txt"):
        #        print item+"/atmos_data.txt does not exist.  Skipping to next file."
        return False
    return True
