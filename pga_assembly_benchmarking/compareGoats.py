#!/usr/bin/python
# Script written by Shawn Sullivan and Phase Genomics
#special purpose script to compare RH map contig orderings to Lachesis goat orderings

import sys
from copy import copy, deepcopy

############################
# USAGE AND INPUT CHECKING #
############################

def printUsage():
    print "\ncompareGoats.py: compare RH map contig orderings to Lachesis goat orderings."
    print "compareGoats.py usage:"
    print "\t./compareGoats.py <rh_map_ordering> <Lachesis_ordering_0> <Lachesis_ordering_1> ...\n"
    return

if len(sys.argv) < 3:
    printUsage()
    sys.exit(1)

###########
# CLASSES #
###########
class Contig():
    """The Contig class stores basic information about a contig necessary to
        facilitate comparison between results in an RH map and results in
        Lachesis ordering files. Specifically, it only requires the name and
        orientation of a contig to do this.
        
        Attributes:
            name(str): the name of the contig
            orientation(str): the orientation of the contig:
                +: forward (relative to the original contig)
                -: reverse complement (relative to the original contig)
                ?: unknown
        """
    
    def __init__(self, name, orientation):
        """Initialize a new Contig"""
        self.name = name.strip()
        self.orientation = orientation
    
    def invert_orientation(self):
        """Invert the orientation of the contig.
            If it's +, it becomes -.
            If it's -, it becomes +.
            If it's ?, it stays ?.
            Otherwise, an Exception is raised.
        """
        if self.orientation == "+":
            self.orientation = "-"
        elif self.orientation == "-":
            self.orientation = "+"
        elif self.orientation == "?":
            self.orientation = "?"
        else:
            raise Exception("Attempting to invert unknown orientation: "
                            + str(self.orientation))

    def equals(self, contig):
        """Return true if self is equal to specified contig, else false."""
        if not isinstance(contig, Contig):
            return False
        return self.name == contig.name

    def __eq__(self, contig):
        """Return true if self is equal to specified contig, else false."""
        if not isinstance(contig, Contig):
            return False
        return self.name == contig.name


    def __repr__(self):
        """Return a string representation of a Contig."""
        return self.name

#############
# FUNCTIONS #
#############

def parse_rh_map(rh_map_ordering):
    """Parses an RH map file and returns a dictionary of
            <chr_n> : [contig1, contig2, ..., contign]
        
        Arguments:
            rh_map_ordering(str): the name of the RH map file to parse
            
        Return:
            dict[str:list[Contig]]: dictionary mapping chromosome numbers (as
                strings) to a list of Contig objects on that chromosome, based
                on parsing the contents of the RH map file.
        
    """
    #we are going to store results in a dictionary
    results = {}
    with open(rh_map_ordering) as f:
        for line in f:
            #line format: chr_n	contig_name	min_loc	max_loc	orientation
            line_split = line.split()
            
            #rh map files list "X" as a chromosome number, which we will call
            #"30" because that's what it's called by Lachesis
            chr_n = line_split[0] if line_split[0] != "X" else "30"
            contig_name = line_split[1]
            orientation = line_split[4]
            
            #make a Contig object to put in our dictionary
            contig = Contig(contig_name, orientation)
            
            #if chr_n isn't in our results yet, create an empty list for it
            if chr_n not in results.keys():
                results[chr_n] = []
            results[chr_n].append(contig)

    return results

def parse_ordering_files(lachesis_orderings):
    """Parses a list of lachesis ordering files and returns a dictionary of
            <chr_n> : [contig1, contig2, ..., contign]
        
        Arguments:
            lachesis_orderings(list[str]): a list of the names of the Lachesis
                ordering files to parse
        
        Return:
            dict[str:list[Contig]]: dictionary mapping chromosome numbers (as
                strings) to a list of Contig objects on that chromosome, based
                on parsing the contents of the Lachesis ordering files.
    """
    def parse_ordering_file(lachesis_ordering):
        #Helper funciton which parses an ordering file and returns a list of
        #   [contig1, contig2, ..., contign]
        results = []
        with open(lachesis_ordering) as f:
            for line in f:
                #line format: contig_ID(local)	contig_name	contig_rc	orientation_Q_score	gap_size_after_contig
                line = line.strip()
                #skip comments
                if len(line) == 0 or line[0] == "#":
                    continue
                line_split = line.split()
                
                #make a Contig object to put in our dictionary
                contig_name = line_split[1]
                orientation = "-" if int(line_split[2]) else "+"
                contig = Contig(contig_name, orientation)
                
                results.append(contig)
    
        return results
    
    #we are going to store results in a dictionary
    results = {}
    #go through each ordering file, get its contigs from the helper function,
    #and accumulate the results in our results dictionary
    chr_n = 1
    for lachesis_ordering in lachesis_orderings:
        results[str(chr_n)] = parse_ordering_file(lachesis_ordering)
        chr_n += 1
    
    return results

def get_hits(lachesis_ordering, rh_ordering, verbose=False):
    """Count the number of contigs in a lachesis ordering that are also present
        in an RH map ordering. (Note that this could work on two arbitrary lists
        of contigs, but variables are named based on intended usage).
        
        Arguments:
            lachesis_ordering(list[Contig]): a list of Contigs in an ordering
                from Lachesis
            rh_ordering(list[Contig]): a list of Contigs in an ordering from an
                RH map
            verbose(bool): whether to print out additional information on which
                specific contigs are hits. Useful mainly for debugging because
                it clutters the output. Default: False
        
        Return:
            int: the number of Contigs in lachesis_ordering that are also in
                rh_ordering
    """
    #loop through Contigs in lachesis_ordering. if a Contig is in the
    #rh_ordering too, increment hits. print what happened if verbose
    hits = 0
    for lachesis_contig in lachesis_ordering:
        if lachesis_contig in rh_ordering:
            hits += 1
            if verbose:
                print lachesis_contig, rh_contig, "YES"
        elif verbose:
                print lachesis_contig, rh_contig, "no"
    return hits

def get_all_matches(lachesis_ordering, rh_results):
    """Get the list of all RH orderings which have at least one hit with a
        specified Lachesis ordering.
        
        Arguments:
            lachesis_ordering(list[Contig]): a list of Contigs in an ordering
                from Lachesis
            rh_results(dict[str:list[Contig]]): a dictionary mapping chromosome
                numbers to Contig objects in an RH map ordering of that
                chromosome. This is the same kind of dictionary generated by
                parse_rh_map
        
        Return:
            list[list[Contig]]: the list of RH orderings which have at least 1
                hit with lachesis_ordering
    """
    #loop through each ordering in rh_results, get the number of hits with
    #lachesis_ordering, and if there is at least 1 hit, append the ordering to
    #the list we will return
    matches = []
    for rh_ordering in rh_results.keys():
        if get_hits(lachesis_ordering, rh_results[rh_ordering]) > 0:
            matches.append(rh_ordering)
    return matches

def get_best_match(lachesis_ordering, rh_results):
    """Get RH orderings which has the most hits with a specified Lachesis
        ordering.
        
        Arguments:
            lachesis_ordering(list[Contig]): a list of Contigs in an ordering
                from Lachesis
            rh_results(dict[str:list[Contig]]): a dictionary mapping chromosome
                numbers to Contig objects in an RH map ordering of that
                chromosome. This is the same kind of dictionary generated by
                parse_rh_map
        
        Return:
            tuple[list[Contig], int]: a tuple containing
                0:  the RH ordering with the most hits with lachesis_ordering
                1:  the number of hits that RH ordering had with
                    lachesis_ordering
        """
    #keep track of the best matching RH order we've seen and how many hits it
    #had with lachesis_ordering
    best_match = None
    best_hits = -1
    #loop through all RH orderings, and when we find one with more hits than
    #we've seen so far, update our tracking variables
    for rh_ordering in rh_results.keys():
        hits = get_hits(lachesis_ordering, rh_results[rh_ordering])
        if hits > best_hits:
            best_hits = hits
            best_match = rh_ordering
    return best_match, best_hits

def is_reverse_orientation_better(fixed_contigs, dependent_contigs):
    """Check whether one set of Contigs would match another better if it were
        reversed.
        
        Arguments:
            fixed_contigs(list[Contig]): the list of Contigs to consider as
                fixed - i.e., don't consider reversing this list of Contigs
            dependent_contigs(list[Contig]): the list of Contigs to consider as
                reversable
        Return:
            bool: True if dependent_contigs would match fixed_contigs better if
                it were reversed; else False (including False if forward and
                reverse are equally good)
    """
    def get_ordered_hits(primary_array, secondary_array):
        #helper function which determines how many elements in primary_array
        #appear in the same order in secondary_array, from left to right,
        #allowing some elements in primary_array to be skipped. equivalent to
        #removing all elements in primary_array that are not also in
        #secondary_array, and then computing the number of elements in
        #secondary_array that appear in the same order as primary_array reading
        #left to right
        
        #use a pointer to keep track of the next element we are trying to match
        #in secondary_array
        secondary_ptr = 0
        hits = 0
        
        #loop over elements in primary array
        for elt in primary_array:
            #if our pointer is past the end of secondary_array, we're done
            if secondary_ptr >= len(secondary_array):
                break
            #if the current element in primary_array matches the element in
            #secondary_array at our pointer, that's a hit
            if secondary_array[secondary_ptr] == elt:
                secondary_ptr += 1
                hits += 1
            #if there was no hit, leave secondary_ptr where it is. we'll see if
            #a later element in primary_array hits it
        return hits

    #compute the number of ordered hits both forward and reverse, and return
    #their comparison
    reverse_dependent_contigs = list(reversed(dependent_contigs))
    forward_count = get_ordered_hits(fixed_contigs, dependent_contigs)
    reverse_count = get_ordered_hits(fixed_contigs, reverse_dependent_contigs)

    return reverse_count > forward_count

def match_order_and_orientation(dependent_contigs, reverse):
    """If reverse is true, reverse a list of Contigs and invert the orientation
        of all Contigs in that list
        
        Arguments:
            dependent_contigs(list[Contig]): the list of Contigs to consider as
                reversable
            reverse(bool): whether or not to reverse dependent_contigs
        
        Return:
            list[Contig]: if reverse is True, dependent_contigs reversed and
                with all Contig orientations inverted. if reverse is False,
                this method does nothing to dependent_contigs and simply returns
                it.
    """
    if reverse:
        reverse_dependent_contigs = list(reversed(dependent_contigs))
        #if we're going to return the reversed contigs, we need to invert all
        #the orientations in the Contigs in it
        for contig in reverse_dependent_contigs:
            contig.invert_orientation()
        return reverse_dependent_contigs

    return dependent_contigs

def get_contig_string(contig_name):
    """Return a nicely formatted string from contig_name, with tabs added such
        that it prints nicely in the report this script makes. This method is
        designed to work with strings of less than 20 characters or so.
        
        Arguments:
            contig_name(str): the string to use as a basis for the nicely
                formatted string
        
        Return:
            str: a nicely formatted string based on contig_name, padded with
                tabs so that strings of different lengths passed into this
                function will be displayed in the same visual area.
    """
    if len(contig_name) < 8:
        return contig_name + "\t\t"
    elif len(contig_name) < 12:
        return contig_name + "\t"
    else:
        return contig_name + ""

def get_contig_match_key(contig, matches_dict):
    """Determine the dictionary key whose value is a list containing a Contig.
        Intended for use with a dictionary like rh_results or lachesis_results.
        
        Arguments:
            contig(Contig): a Contig to search for
            matches_dict(dict[str:list[Contig]]): a dictionary to search.
        
        Return:
            str: the key whose value contains contig, or "-1" if no such key
                is found
    """
    for key in matches_dict.keys():
        if contig in matches_dict[key]:
            return key
    return "-1"

def is_contig_in_matches(contig, matches_dict):
    """Determine whether a given Contig is inside the lists in the values of a
        dictionary. Intended for use with a dictionary like rh_results or
        lachesis_results.
        
        Arguments:
            contig(Contig): a Contig to search for
            matches_dict(dict[str:list[Contig]]): a dictionary to search.
        
        Return:
            bool: True if Contig can be found inside one of the lists in the
                dictionary's keys; else False
        """
    return get_contig_match_key(contig, matches_dict) != "-1"

def get_contig_match_value(contig, matches_dict):
    """Get the value containing a given Contig from a dictionary. Intended for
        use with a dictionary like rh_results or lachesis_results.
        
        Arguments:
            contig(Contig): a Contig to search for
            matches_dict(dict[str:list[Contig]]): a dictionary to search.
        
        Return:
        list[Contig]: the list of Contigs from matches_dict's values containing
            contig, or None if no such value exists
        """
    key = get_contig_match_key(contig, matches_dict)
    if key != "-1":
        return matches_dict[key]
    else:
        return None

def get_rh_hit_dict(lachesis_ordering, rh_results):
    """Get a dictionary containing information about the hits between a Lachesis
        ordering and a set of RH results. Specifically, get a dictionary mapping
        all chromosomes in rh_results to the list of Contigs found both in the
        values of rh_results for that chromosome and in lachesis_ordering.
        
        Arguments:
            lachesis_ordering(list[Contig]): a list of Contigs in an ordering
                from Lachesis
            rh_results(dict[str:list[Contig]]): a dictionary mapping chromosome
                numbers to Contig objects in an RH map ordering of that
                chromosome. This is the same kind of dictionary generated by
                parse_rh_map
        
        Return:
            dict[str:list[Contig]]: a dictionary mapping keys from rh_results
                (typically chromosome numbers inside strings) to a list of
                contigs that were found in lachesis_ordering and the value in
                rh_results corresponding to each key
    """
    #go through each key (==chromosome) in rh_results. then go through each
    #contig in lachesis_ordering. if the values for the key contain the contig,
    #store this as a hit for that chromosome by storing it in a new dictionary.
    hit_dict = {}
    for key in rh_results.keys():
        hit_candidates = rh_results[key]
        for lachesis_contig in lachesis_ordering:
            if lachesis_contig in hit_candidates:
                if key not in hit_dict.keys():
                    hit_dict[key] = []
                hit_dict[key].append(hit_candidates[hit_candidates.index(lachesis_contig)]) #because it will actually have a different contig object, with the same name but different potentially orientation
    return hit_dict

def print_ordering_comparison_details(rh_results, lachesis_results):
    """Print additional detailed information comparing rh_results to
        lachesis_results. Primarily, prints whether the contigs in each set of
        results are on the same chromosome or different chromosomes in those
        results. Returns nothing.
        
        Arguments:
        rh_results(dict[str:list[Contig]]): a dictionary mapping chromosome
            numbers to Contig objects in an RH map ordering of that
            chromosome. This is the same kind of dictionary generated by
            parse_rh_map
        lachesis_results(dict[str:list[Contig]]): a dictionary mapping
            chromosome numbers to Contig objects in a set of Lachesis orderings.
            This is the same kind of dictionary generated by
            parse_ordering_files
    """
    if len(rh_results.keys()) != len(lachesis_results.keys()):
        print "WARNING: different number of chromosomes detected! RH: " + str(len(rh_results.keys())) + ", LACH: " + str(len(lachesis_results.keys()))

    for i in range(0, max(len(rh_results.keys()), len(lachesis_results.keys()))):
        hit_rh = 1 if i < len(rh_results.keys()) else 0
        hit_lachesis = 1 if i < len(lachesis_results.keys()) else 0
        chr_n = i+1
        if hit_rh and hit_lachesis:
            print "Hit chr " + str(chr_n) + " for both RH and Lachesis!"
            print "Ordering:"
            rh_contigs = rh_results[str(chr_n)]
            lachesis_contigs = lachesis_results[str(chr_n)]
            for j in range(0, max(len(rh_contigs), len(lachesis_contigs))):
                print_str = str(j) + "\t"
                print_str = print_str + str(rh_contigs[j].name) + "\t" if j < len(rh_contigs) else print_str + "\t\t"
                print_str = print_str + str(lachesis_contigs[j].name) if j < len(lachesis_contigs) else print_str
                print print_str
        elif hit_rh:
            print "Hit chr " + str(chr_n) + " for RH only!"
            print "Ordering:"
            rh_contigs = rh_results[str(chr_n)]
            for j in range(0, len(rh_contigs)):
                print_str = str(j) + "\t" + str(rh_contigs[j].name)
                print print_str
        elif hit_lachesis:
            print "Hit chr " + str(chr_n) + " for Lachesis only!"
            print "Ordering:"
            lachesis_contigs = lachesis_results[str(chr_n)]
            for j in range(0, len(lachesis_contigs)):
                print_str = str(j) + "\t" + str(lachesis_contigs[j].name)
                print print_str
        else:
            print "No hits for chr " + str(chr_n) + "!"

def print_best_match_comparison_details(rh_results, lachesis_results):
    """Print additional detailed information comparing rh_results to
        lachesis_results. Primarily, prints information about the best matches
        for each contig in each set of results, and whether or not the two sets
        of results agree. Returns nothing.
        
        Arguments:
            rh_results(dict[str:list[Contig]]): a dictionary mapping chromosome
                numbers to Contig objects in an RH map ordering of that
                chromosome. This is the same kind of dictionary generated by
                parse_rh_map
            lachesis_results(dict[str:list[Contig]]): a dictionary mapping
                chromosome numbers to Contig objects in a set of Lachesis orderings.
                This is the same kind of dictionary generated by
                parse_ordering_files
        """
    print "\n===============\nBest match comparison"

    for i in range(1, len(lachesis_results.keys())+1):
        print "LACHESIS CHR" + str(i-1)
        print "Contig\t\t\tLp\tRHp\tLo\tRHo"
        lachesis_ordering = lachesis_results[str(i)]
        best_match_name = get_best_match(lachesis_ordering, rh_results)[0]
        reverse = is_reverse_orientation_better(lachesis_ordering, rh_results[best_match_name])
        best_results = match_order_and_orientation(rh_results[best_match_name], reverse)

        lachesis_ptr = 0
        rh_ptr = 0
        while True:
            lachesis_contig = lachesis_ordering[lachesis_ptr] if lachesis_ptr < len(lachesis_ordering) else None
            rh_contig = best_results[rh_ptr] if rh_ptr < len(best_results) else None
            rh_ptr_display = str(rh_ptr)
            if reverse:
                rh_ptr_display = str(len(best_results) - rh_ptr)
            #check if we're done
            if lachesis_contig is None and rh_contig is None:
                break
            #if the contigs match
            elif lachesis_contig == rh_contig:
                lachesis_ptr += 1
                rh_ptr += 1
                print get_contig_string(lachesis_contig.name) + "\t" + str(lachesis_ptr) + "\t" + rh_ptr_display + "\t" + str(lachesis_contig.orientation) + "\t" + str(rh_contig.orientation)
            #if we've already run past the end of the rh map contigs
            elif rh_contig is None:
                lachesis_ptr += 1
                print get_contig_string(lachesis_contig.name) + "\t" + str(lachesis_ptr) + "\t.\t" + str(lachesis_contig.orientation) + "\t."
            #if we've already run past the end of the lachesis contigs
            elif lachesis_contig is None:
                rh_ptr += 1
                print get_contig_string(rh_contig.name) + "\t.\t" + rh_ptr_display + "\t.\t" + str(rh_contig.orientation)
            #if the lachesis contig is later on in the rh contigs
            elif lachesis_contig in best_results[rh_ptr+1:]:
                rh_ptr += 1
                print get_contig_string(rh_contig.name) + "\t.\t" + rh_ptr_display + "\t.\t" + str(rh_contig.orientation)
            #if the rh contig is later on in the lachesis contigs
            elif rh_contig in lachesis_ordering[lachesis_ptr+1:]:
                lachesis_ptr += 1
                print get_contig_string(lachesis_contig.name) + "\t" + str(lachesis_ptr) + "\t.\t" + str(lachesis_contig.orientation) + "\t."
            #if none of the above, print the lachesis contig if available, otherwise the rh contig
            elif lachesis_contig is not None:
                lachesis_ptr += 1
                print get_contig_string(lachesis_contig.name) + "\t" + str(lachesis_ptr) + "\t.\t" + str(lachesis_contig.orientation) + "\t."
            elif rh_contig is not None:
                rh_ptr += 1
                print get_contig_string(rh_contig.name) + "\t.\t" + rh_ptr_display + "\t.\t" + str(rh_contig.orientation)
            #should never happen
            else:
                raise Exception("ERROR: can't handle contigs. Lachesis=" + str(lachesis_contig) + ", RH=" + str(rh_contig))
        print

    print "\n===============\nBest match comparison"

    for i in range(1, len(lachesis_results.keys())+1):
        print "LACHESIS CHR" + str(i-1)
        print "Contig\t\t\tLp\tRHp\tLo\tRHo"
        lachesis_ordering = lachesis_results[str(i)]
        best_match_name = get_best_match(lachesis_ordering, rh_results)[0]
        best_results = copy(rh_results[best_match_name])
        leftovers = copy(best_results)

        for j in range(0, len(lachesis_ordering)):
            lachesis_contig = lachesis_ordering[j]
            if lachesis_contig in best_results:
                rh_contig = best_results[best_results.index(lachesis_contig)]   #because it will actually have a different contig object, with the same name but different potentially orientation
                print get_contig_string(lachesis_contig.name) + "\t" + str(j) + "\t" + str(best_results.index(rh_contig)) + "\t" + str(lachesis_contig.orientation) + "\t" + str(rh_contig.orientation)
                leftovers.remove(rh_contig)
            else:
                print get_contig_string(lachesis_contig.name) + "\t" + str(j) + "\t.\t" + str(lachesis_contig.orientation) + "\t."
        #print the leftovers
        for rh_contig in leftovers:
            rh_contig_display = rh_contig.name
            if rh_contig in lachesis_ordering:
                rh_contig_display = "*" + rh_contig_display
                #each time through when we hit a contig present multiple times (that is an assumption at this point in the code - if something is in leftovers and
                #lachesis_ordering still, it must have appeared more than once), we need to clear out one occurance from best_results so we get the right index
                #this depends on best_results being a copy of the actual results
                best_results[best_results.index(rh_contig)] = None
            print get_contig_string(rh_contig_display) + "\t.\t" + str(best_results.index(rh_contig)) + "\t.\t" + str(rh_contig.orientation)
        print



###############
# MAIN METHOD #
###############

#enabling verbose causes print_ordering_comparison_details and
#print_best_match_comparison_details to be called
verbose = False

#parse command line options, very fragily
rh_map_ordering = sys.argv[1]
lachesis_orderings = sys.argv[2:]

rh_results = parse_rh_map(rh_map_ordering)
lachesis_results = parse_ordering_files(lachesis_orderings)

if verbose:
    print_ordering_comparison_details(rh_results, lachesis_results)


#print stats for how Lachesis results compare to RH map results. this comparison
#is one directional. we will find out if Lachesis said something that the RH
#map differed with, but we won't find out if Lachesis *didn't* say something the
#RH map said - we do that later.
print "Lachesis -> RH map hits"

for i in range(1, len(rh_results.keys())+1):
    if i != 30:
        print "\t" + str(i),
    else:
        print "\tX\tN"
print
unique_maps = []
non_unique_maps = []
all_ctgs_hits = 0
all_best_ctgs_hits = 0
all_ctgs = 0
for i in range(1, len(lachesis_results.keys())+1):
    hitters = 0
    print_str = str(i) + "\t"
    total_hits = 0
    for j in range(1, len(rh_results.keys())+1):
        hits = get_hits(lachesis_results[str(i)], rh_results[str(j)])
        total_hits += hits
        if hits > 0:
            hit_percent = round((hits/float(len(lachesis_results[str(i)]))) * 100.0, 1)
            if hit_percent < 10:
                print_str += str(hit_percent) + "%"
            elif hit_percent == 100:
                print_str += str(int(hit_percent)) + "%"
            else:
                print_str += str(int(hit_percent)) + "%\t"
            hitters += 1
        else:
            print_str += ".\t"
    all_ctgs_hits += total_hits
    all_ctgs += len(lachesis_results[str(i)])
    best_match, best_hits = get_best_match(lachesis_results[str(i)], rh_results)
    all_best_ctgs_hits += best_hits
    miss_percent = round(((len(lachesis_results[str(i)])-total_hits)/float(len(lachesis_results[str(i)]))) * 100.0, 1)
    if miss_percent == 0:
        print_str += "."
    else:
        print_str += str(miss_percent) + "%"
    print print_str
    if hitters == 1:
        unique_maps.append(i)
    else:
        non_unique_maps.append(i)

print "\nUnique hitters: " + str(round(len(unique_maps)/float(len(lachesis_results.keys())) * 100.0, 1)) + "%"
print unique_maps
print "Non-unique hitters: " +  str(round(len(non_unique_maps)/float(len(lachesis_results.keys())) * 100.0, 1)) + "%"
print non_unique_maps
print "Total hits: " + str(all_ctgs_hits) + " of " + str(all_ctgs) + ", " + str(round(all_ctgs_hits/float(all_ctgs) * 100.0, 1)) + "%"
print "Best hits: " + str(all_best_ctgs_hits) + " of " + str(all_ctgs) + ", " + str(round(all_best_ctgs_hits/float(all_ctgs) * 100.0, 1)) + "%"

#print stats for how RH map results compare to Lachesis results. this comparison
#is one directional in a way complementary to the previous printing.
print "\n===============\nRH map -> Lachesis hits"

for i in range(1, len(lachesis_results.keys())+1):
    print "\t" + str(i),
print "\tN"
unique_maps = []
non_unique_maps = []
full_hitters = []
non_full_hitters = []
unique_full_hitters = []
non_unique_full_hitters = []
all_ctgs_hits = 0
all_best_ctgs_hits = 0
all_ctgs = 0

for i in range(1, len(rh_results.keys())+1):
    hitters = 0
    print_str = str(i) + "\t"
    if i == 30:
        print_str = "X\t"
    total_hits = 0
    for j in range(1, len(lachesis_results.keys())+1):
        hits = get_hits(rh_results[str(i)], lachesis_results[str(j)])
        total_hits += hits
        if hits > 0:
            hit_percent = round((hits/float(len(rh_results[str(i)]))) * 100.0, 1)
            if hit_percent < 10:
                print_str += str(hit_percent) + "%"
            elif hit_percent == 100:
                print_str += str(int(hit_percent)) + "%"
            else:
                print_str += str(int(hit_percent)) + "%\t"
            hitters += 1
        else:
            print_str += ".\t"
    all_ctgs_hits += total_hits
    all_ctgs += len(rh_results[str(i)])
    best_match, best_hits = get_best_match(lachesis_results[str(i)], rh_results)
    all_best_ctgs_hits += best_hits
    miss_percent = round(((len(rh_results[str(i)])-total_hits)/float(len(rh_results[str(i)]))) * 100.0, 1)
    if miss_percent == 0:
        print_str += "."
    else:
        print_str += str(miss_percent) + "%"
    print print_str
    if hitters == 1:
        unique_maps.append(i)
    else:
        non_unique_maps.append(i)
    if total_hits == len(rh_results[str(i)]):
        full_hitters.append(i)
    else:
        non_full_hitters.append(i)
    if total_hits == len(rh_results[str(i)]) and hitters == 1:
        unique_full_hitters.append(i)
    else:
        non_unique_full_hitters.append(i)

def replace_30_with_X(array):
    if 30 in array:
        array[array.index(30)] = "X"

replace_30_with_X(unique_maps)
replace_30_with_X(non_unique_maps)
replace_30_with_X(full_hitters)
replace_30_with_X(non_full_hitters)
replace_30_with_X(unique_full_hitters)
replace_30_with_X(non_unique_full_hitters)

print "\nUnique hitters: " + str(round(len(unique_maps)/float(len(rh_results.keys())) * 100.0, 1)) + "%"
print unique_maps
print "Non-unique hitters: " +  str(round(len(non_unique_maps)/float(len(rh_results.keys())) * 100.0, 1)) + "%"
print non_unique_maps
print "Full hitters: " + str(round(len(full_hitters)/float(len(rh_results.keys())) * 100.0, 1)) + "%"
print full_hitters
print "Non-full hitters: " + str(round(len(non_full_hitters)/float(len(rh_results.keys())) * 100.0, 1)) + "%"
print non_full_hitters
print "Unique-full hitters: " + str(round(len(unique_full_hitters)/float(len(rh_results.keys())) * 100.0, 1)) + "%"
print full_hitters
print "Non-unique-full hitters: " + str(round(len(non_unique_full_hitters)/float(len(rh_results.keys())) * 100.0, 1)) + "%"
print non_full_hitters
print "Total hits: " + str(all_ctgs_hits) + " of " + str(all_ctgs) + ", " + str(round(all_ctgs_hits/float(all_ctgs) * 100.0, 1)) + "%"
print "Best hits: " + str(all_best_ctgs_hits) + " of " + str(all_ctgs) + ", " + str(round(all_best_ctgs_hits/float(all_ctgs) * 100.0, 1)) + "%"


#next, instead of just getting the best match for each contig, get all the
#matches
#then print out the Lachesis contigs and where they can be found in the RH maps
#show places where the Lachesis contig isn't placed anywhere by the RH maps
#also show where the RH maps have things Lachesis doesn't if that happens
#show whether the orientations agree or not; if we "flip" the RH map to line
#things up, we may need to flip the orientation
#example:
#LACHESIS CHR0
#contig1    rh_chr.position orientation_concordance
#contig2    rh_chr.position orientation_concordance
#contig3                    orientation
#contig4                    orientation
#contig5    rh_chr.position orientation_concordance
#           rh_chr.position orientation
#contig6    rh_chr.position orientation_concordance
#etc

if verbose:
    print_best_match_comparison_details(rh_results, lachesis_results)

print "\n===============\nAll Lachesis->RH map matches comparison"

all_up_hits = 0
all_up_misses = 0
all_up_false_misses = 0

for i in range(1, len(lachesis_results.keys())+1):
    print "LACHESIS CHR" + str(i-1)
    print "Contig\t\t\tLp\tRHc\tRHp\tLo\tRHo"
    lachesis_ordering = lachesis_results[str(i)]
    matches = get_all_matches(lachesis_ordering, rh_results)
    match_results = {}
    for match in matches:
        match_results[match] = copy(rh_results[match])

    match_leftovers = get_rh_hit_dict(lachesis_ordering, rh_results)
    hit_contigs = []

    total_hits = 0
    total_misses = 0
    false_misses = 0
    for j in range(0, len(lachesis_ordering)):
        lachesis_contig = lachesis_ordering[j]
        if is_contig_in_matches(lachesis_contig, match_results):
            total_hits += 1
            hit_contigs.append(lachesis_contig)
            rh_contig_key = get_contig_match_key(lachesis_contig, match_results)
            rh_contig_value = get_contig_match_value(lachesis_contig,
                                                     match_results)
            #because it will actually have a different contig object, with the
            #same name but different potentially orientation
            rh_contig = rh_contig_value[rh_contig_value.index(lachesis_contig)]
            print get_contig_string(lachesis_contig.name) + "\t" + str(j) + "\t"\
                    + str(rh_contig_key) + "\t"\
                    + str(rh_contig_value.index(rh_contig)) + "\t"\
                    + str(lachesis_contig.orientation) + "\t"\
                    + str(rh_contig.orientation)
            match_leftovers[rh_contig_key].remove(rh_contig)
        else:
            print get_contig_string(lachesis_contig.name) + "\t" + str(j)\
                + "\t.\t.\t" + str(lachesis_contig.orientation) + "\t."

    #print the leftovers: things we didn't find a match for
    for key in match_leftovers:
        leftovers = match_leftovers[key]
        for rh_contig in leftovers:
            rh_contig_display = rh_contig.name
            if rh_contig in hit_contigs:
                false_misses += 1
                rh_contig_display = "*" + rh_contig_display
                #each time through when we hit a contig present multiple times
                #(that is an assumption at this point in the code - if something
                #is in leftovers and lachesis_ordering still, it must have
                #appeared more than once), we need to clear out one occurance
                #from best_results so we get the right index this depends on
                #best_results being a copy of the actual results
                leftovers[leftovers.index(rh_contig)] = None
            else:
                total_misses += 1
            print get_contig_string(rh_contig_display) + "\t.\t" + key + "\t"\
                + str(rh_results[key].index(rh_contig)) + "\t.\t"\
                + str(rh_contig.orientation)
    all_up_hits += total_hits
    all_up_misses += total_misses
    all_up_false_misses += false_misses
    print "Summary: " + str(total_hits) + " possible RH matches hit, "\
        + str(round(float(total_hits)/(total_hits + total_misses + false_misses) * 100.0, 1)) + "%."
    if total_misses + false_misses != 0:
        print "Of the " + str(total_misses + false_misses) + " misses, "\
            + str(false_misses) + " were multiple-mapping RH matches. Total hits plus false matches is "\
            + str(total_hits + false_misses) + ", "\
            + str(round(float(total_hits+false_misses)/(total_hits + total_misses + false_misses) * 100.0, 1)) + "%"
    print

#final summary
print "ALL UP SUMMARY"
print "Summary: " + str(all_up_hits) + " possible RH matches hit, "\
    + str(round(float(all_up_hits)/(all_up_hits + all_up_misses + all_up_false_misses) * 100.0, 1)) + "%."
if all_up_misses + all_up_false_misses != 0:
    print "Of the " + str(all_up_misses + all_up_false_misses) + " misses, "\
        + str(all_up_false_misses) + " were multiple-mapping RH matches. Total hits plus false matches is "\
        + str(all_up_hits + all_up_false_misses) + ", "\
        + str(round(float(all_up_hits+all_up_false_misses)/(all_up_hits + all_up_misses + all_up_false_misses) * 100.0, 1)) + "%"


