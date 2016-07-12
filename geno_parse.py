#!/usr/bin/python3
from collections import defaultdict
import bisect
import argparse

# if any portion of any range in to_remove is in
def read_gaps(gap_file, ranges_to_remove):
    """Reads ranges from a gapfile

    Adds each range to the ranges_to_remove parameter. Can be used successively
    to read several files.

    Args:
        gap_file: A file containing 'gap' data. One white-space separated tuple
                of the form 'chromosome_name gap_begin gap_end' is expected
                per line. Tokens after the first 3 are ignored.
                Improperly formatted lines are ignored.
        ranges_to_remove: A defaultdict in which the gaps will be stored:
                ranges_to_remove[chromosome_name]=[(gap_begin, gap_end),...]
    """
    # make dict of items to remove {chr: [(start, stop)]}
    with open(gap_file) as f:
        for line in f:
            line = line.split()[:3]
            try:
                bisect.insort(ranges_to_remove[line[0]],
                        (int(line[1]), int(line[2])))
            except ValueError:
                print ('-->Line not added to ranges_to_remove: {}'
                        .format(' '.join(line)))


def merge_ranges(ranges):
    """Merges adjacent ranges in the ranges_to_remove dictionary.

    Note: Ranges are INCLUSIVE of begin- AND end-points.

    Args:
        ranges: A defultdict containing lists of ranges for a set of
                chromosomes.
    """
    for geno in ranges:
        i = 0
        while i < len(ranges[geno]) - 1:
            end = i
            while (end < len(ranges[geno]) - 1 and
                    ranges[geno][end][1] >= ranges[geno][end+1][0] - 1):
                end += 1
            ranges[geno][i] = (ranges[geno][i][0], ranges[geno][end][1])
            # remove merged ranges
            for j in range(i, end):
                del ranges[geno][i+1]
            i += 1



        # for i, t in enumerate(ranges[geno]):
        #     if i < len(ranges[geno]) and t[1] >= ranges[geno][i+1][0] - 1:
        #         ranges[geno][i] = (t[0], ranges[geno][i+1][1])
        #         del ranges[geno][i+1]


def filter_ranges(remaining, ranges, max_overlap=0):
    """Divides a list of input lines into two: a list of those that were not
    removed due to an overlap with a range from a gap_file, and those that were.

    'Track' heading lines are duplicated when component lines are present in the
    list of removed lines.

    Args:
        remaining: A list of lines from the input bed-file. It is expected that
                the first 3 white-space separated tokens in each data line will
                be chromosome_name, range_begin, and range_end.
        ranges: A defaultdict of ranges with which a matching overlap will cause
                removal of a data line from the 'remaining' list.
        max_overlap: The maximum length intersection allowed before a data line
                is removed.

    Returns:
        removed: A list of removed lines including duplicated 'track' header.
    """
    removed_lines = []
    removed_idxs = []
    track = ''

    for i, line in enumerate(remaining):
        line = line.strip()
        # save track lines for necessary entries in removed list
        if line.startswith('track'):
            track = line
            continue
        # check all other lines for overlap with remove ranges
        else:
            geno, start, stop = line.split()[:3]
            start, stop = int(start), int(stop)
            overlap = range(0)
            # check this line against each gap
            for t in ranges[geno]:
                overlap = range(max(start, t[0]), min(stop, t[1])+1)
                if len(overlap) > max_overlap or stop < t[0]:
                    break
            if len(overlap) > max_overlap:
                # add track line and removed line to removed_lines list
                if track:
                    removed_lines.append(track)
                    track = ''
                removed_lines.append(line)
                # add this index to a list of lines to be removed
                removed_idxs.append(i)
    # delete removed lines from remaining list
    for i in range(len(removed_idxs)-1, -1, -1):
        del remaining[removed_idxs[i]]
    return removed_lines

def main():
    # argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', required=True,
            help=('File from which lines with overlapping ranges are to be'
                  ' removed'))
    parser.add_argument('-g', '--gap_files', nargs='+', required=True,
            help=('Space separated list of files listing gaps to remove from'
                  ' infile'))
    parser.add_argument('-o', '--outfile', nargs='?', default='filtered.bed',
            help=('Optional output file name for retained lines; defalut:'
                  'filtered.bed'))
    parser.add_argument('-r', '--removed', nargs='?',
            help=('Enable output of removed lines to separate file, by'
                  ' specifying a second output file name.'))
    parser.add_argument('-m', '--max_overlap', nargs='?', default=0,
            help=('Maximum allowed intersection size of input line and gap'
                  ' ranges.'))
    args = parser.parse_args()

    # read gap files into dict of sorted lists
    ranges_to_remove = defaultdict(list)
    for gap_file in args.gap_files:
        print ('Reading {}...'.format(gap_file))
        read_gaps(gap_file, ranges_to_remove)
        print ('Done reading {}.'.format(gap_file))

    # merge any overlapping ranges in dict of ranges to remove
    print ('Merging ranges to remove...')
    merge_ranges(ranges_to_remove)
    print ('Done merging.')

    # read input file
    with open(args.infile) as i:
        remaining = [line for line in i]
    # filter input to separate overlapping and non-overlapping lines
    removed = filter_ranges(remaining, ranges_to_remove, args.max_overlap)
    # write lines with no overlap
    with open(args.outfile, 'w') as o:
        for line in remaining:
            print(line, file=o, end='')
    # write removed lines to separate file
    if args.removed:
        with open(args.removed, 'w') as r:
            for l in removed:
                print(l, file=r)

if __name__ == '__main__':
    main()
