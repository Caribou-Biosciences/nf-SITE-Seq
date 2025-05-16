#########################################################################
#########################################################################
##            Classes for manipulating genomic coordinates             ##
#########################################################################
#########################################################################


import re
from abc import abstractmethod


class Coordinate(object):
    __slots__ = ["chrom", "start", "stop", "is_reverse_strand"]

    _pattern = re.compile(
        r"^(?P<chrom>[A-Za-z0-9-_.]+):(?P<start>\d+)(?:[-_](?P<stop>\d+))?,?(?P<strand>[+-])?$"
    )

    @abstractmethod
    def __init__(self, chrom, start=None, stop=None, is_reverse_strand=None):
        """
        Chromosome coordinate abstract class. Onedinate and Zerodinate will extend this abstract class to 1- and 0-
        based chromosome coordinate. The coordinate start will always <= stop
        chrom must be provided. Coordinates can be provided through each variable or by sending to chrom altogether.
        If is_reverse_strand is not provided, then it is assumed to be False.
        Valid input:
            Coordinate('chr1', 123, 456, False)
            Coordinate('chr1:123-456,+')
            Coordinate('chr1:123,+')
            Coordinate('chr1:123')
        Not valid:
            Coordinate('chr1', is_reverse_strand=False) : missing start_position
            Coordinate('chr1:123', is_reverse_strand=False) : using string and separate params to pass in coordinates
            simultaneously is not allowed.
        Args:
            chrom (str): chromosome
            start (int or None): start_localhost_server position 1-based, should always <=stop
            stop (int or None): stop position 1-based, should always >=start
            is_reverse_strand (bool or str or None): is on reverse strand? or '+' for positive strand, '-' for reverse
            strand
        """
        if (
            (chrom is not None)
            and (start is None)
            and (stop is None)
            and (is_reverse_strand is None)
        ):
            chrom, start, stop, is_reverse_strand = self.parse_string(chrom)

        if (chrom is not None) and (start is not None):
            _chr = str(chrom)
            _start = int(start)
            if stop is None:
                _stop = _start
            elif isinstance(stop, bool):
                raise Exception("Stop cannot be boolean")
            else:
                _stop = int(stop)

            if isinstance(is_reverse_strand, bool):
                _is_reverse_strand = is_reverse_strand
            elif isinstance(is_reverse_strand, str) and is_reverse_strand == "-":
                _is_reverse_strand = True
            else:
                _is_reverse_strand = False
        else:
            raise Exception("Input cannot be parsed")

        if _start > _stop:
            if _is_reverse_strand is False:
                raise Exception(
                    "Invalid start position (start cannot be larger than the stop)"
                )
            else:
                _start, _stop = _stop, _start
                _is_reverse_strand = True
        elif _start < 0:
            print("WARNING: negative coordinates found, which is unusual")

        super().__setattr__("chrom", _chr)
        super().__setattr__("start", _start)
        super().__setattr__("stop", _stop)
        super().__setattr__("is_reverse_strand", _is_reverse_strand)

    @property
    def chr(self):
        return self.chrom

    @classmethod
    def parse_string(cls, coordinate):
        """
        Process a string into Onedinate variables
        Args:
            coordinate (str): a string representing a Onedinate
        Returns:
            (tuple)
            The extracted string value or None of chrom, start, stop, strand
        """

        matched = cls._pattern.match(coordinate)
        if matched:
            chrom, start, stop, strand = matched.groups()
        else:
            raise Exception("Chromosome coordinate cannot be extracted")
        return str(chrom), start, stop, strand

    def switch_strand(self):
        """
        Change this Onedinate to the other stand (toggle is_reverse_strand)
        Returns:
            (Coordinate): a new Coordinate object with a different strand
        """

        return type(self)(self.chrom, self.start, self.stop, not self.is_reverse_strand)

    def extend_5prime(self, length):
        """
        Extend the 5 prime end of the length nucleotides
        Args:
            length (int): length
        Returns:
            (Coordinate): a new Coordinate object with a extended 5 prime length
        """

        if self.is_reverse_strand:
            return type(self)(self.chrom, self.start, self.stop + length, True)
        else:
            return type(self)(self.chrom, self.start - length, self.stop, False)

    def extend_3prime(self, length):
        """
        Extend the 3 prime end of the length nucleotides.
        Args:
            length (int): length
        Returns:
            (Coordinate): a new Coordinate object with a extended 3 prime length
        """
        if self.is_reverse_strand:
            return type(self)(self.chrom, self.start - length, self.stop, True)
        else:
            return type(self)(self.chrom, self.start, self.stop + length, False)

    def clip_5prime(self, length):
        """
        Clip the 5 prime end of this coordinate by the given length.
        Args:
            length (int): length to clip.
        Returns:
            (Coordinate): a new Coordinate object with 5 prime end clipped.
        """
        if length >= self.length:
            raise Exception(
                f"Cannot clip length {length} when total length is {self.length}"
            )
        elif self.is_reverse_strand:
            return type(self)(self.chrom, self.start, self.stop - length, True)
        else:
            return type(self)(self.chrom, self.start + length, self.stop, False)

    def clip_3prime(self, length):
        """
        Clip the 3 prime end of this coordinate by the given length.
        Args:
            length (int): length to clip.
        Returns:
            (Coordinate): a new Coordinate object with 3 prime end clipped.
        """
        if length >= self.length:
            raise Exception(
                f"Cannot clip length {length} when total length is {self.length}"
            )
        elif self.is_reverse_strand:
            return type(self)(self.chrom, self.start + length, self.stop, True)
        else:
            return type(self)(self.chrom, self.start, self.stop - length, False)

    def intersection(self, coord, ignore_strand=False):
        """
        Return the coordinate of the intersection between this coordinate and
        coord if it exists, None otherwise. If ignore_strand is True, coordinates
        must be on the same strand to intersect. Otherwise, the intersection will
        be returned on the same strand as self (if it exists).
        Args:
            coord (Coordinate): other coordinate to find intersection with.
            Must be the same type.
            ignore_strand (bool, default False): whether the coordinates must be
            on the same strand to intersect. If False, the intersection on the
            same strand as self will be returned (if it exists).
        Returns:
            (Coordinate or None): a new Coordinate object representing the
            intersection of the two coordinates if it exists, None otherwise.
        """
        self_type = type(self)
        if not isinstance(coord, self_type):
            raise Exception(
                f"Cannot find intersection between self of type {self_type} and "
                f"coord of type {type(coord)}"
            )

        if ignore_strand and self.is_reverse_strand != coord.is_reverse_strand:
            coord = coord.switch_strand()

        if not self.is_overlapping_with(coord):
            return None

        return self_type(
            self.chrom,
            max(self.start, coord.start),
            min(self.stop, coord.stop),
            self.is_reverse_strand,
        )

    @classmethod
    def parse_cut_site_string(cls, cut_site_str):
        """
        Parse a formatted cut site string and return a list of cut locations.
        Cut site strings are formatted as such: 'chr1:123,+;chr1:123,-;'
        Args:
            cut_site_str (str):
                A formatted cut site string.
        Returns:
            A list of cut site coordinates.
        """
        cut_locations = filter(None, cut_site_str.split(";"))
        cut_locations = sorted(set([site.split(",")[0] for site in cut_locations]))

        return list(cls(location) for location in cut_locations)

    def serialize(self):
        """
        Return a dictionary that contains the attributes of this Coordinate
        """
        return {
            "chrom": self.chrom,
            "start": self.start,
            "stop": self.stop,
            "is_reverse_strand": self.is_reverse_strand,
        }

    @classmethod
    def deserialize(cls, sdata):
        """
        Return a Coordinate based on the input dictionary produced by the serialize method.
        """
        return cls(
            sdata["chrom"], sdata["start"], sdata["stop"], sdata["is_reverse_strand"]
        )

    @staticmethod
    def combine_coordinates(coordinates):
        """
        This method takes in a list of coordinates and returns a coordinate representing the full range
        of the input coordinates. All coordinates must be of the same type and on the same chromosome;
        if all coordinates are on the same strand the returned coordinate will be on that strand, and
        will otherwise default to the positive strand.
        """
        if len(coordinates) < 2:
            raise Exception("At least two coordinates must be given to combine.")

        same_strand = True
        for coordinate_1, coordinate_2 in zip(coordinates[:-1], coordinates[1:]):
            if type(coordinate_1) is type(coordinate_2):
                raise Exception("Given coordinates are not of the same type.")
            elif coordinate_1.chrom != coordinate_2.chrom:
                raise Exception("Given coordinates are not on the same chromosome.")
            elif coordinate_1.is_reverse_strand != coordinate_2.is_reverse_strand:
                same_strand = False

        all_positions = [
            pos for chroord in coordinates for pos in [chroord.start, chroord.stop]
        ]
        min_position, max_position = min(all_positions), max(all_positions)
        chrom, chroord_type = coordinates[0].chrom, type(coordinates[0])
        if same_strand:
            is_reverse_strand = coordinates[0].is_reverse_strand
        else:
            is_reverse_strand = False

        return chroord_type(
            chrom=chrom,
            start=min_position,
            stop=max_position,
            is_reverse_strand=is_reverse_strand,
        )

    @property
    def boundary_position_5prime(self):
        if self.is_reverse_strand:
            return self.stop
        else:
            return self.start

    @property
    def boundary_position_3prime(self):
        if self.is_reverse_strand:
            return self.start
        else:
            return self.stop

    @property
    @abstractmethod
    def length(self):
        pass

    @property
    @abstractmethod
    def string_value(self):
        pass

    @property
    @abstractmethod
    def string_value_unstranded(self):
        pass

    def __str__(self):
        return self.string_value

    def __repr__(self):
        return self.string_value

    # followings will make dictionary key work and will make set work and will make == and != work
    def __hash__(self):
        return hash(self.string_value)

    def __eq__(self, other):
        return self.string_value == other.string_value

    def __ne__(self, other):
        return self.string_value != other.string_value


class Onedinate(Coordinate):
    def __init__(self, chrom, start=None, stop=None, is_reverse_strand=None):
        super(Onedinate, self).__init__(chrom, start, stop, is_reverse_strand)

    def is_overlapping_with(self, target, ignore_strand=False):
        """
        Check if the current Onedinate is overlapping with the target one.
        The strand does affect the result. If ignore_strand is True,
        coordinates must be on the same strand to overlap.
        Args:
            target (Onedinate): Another Onedinate to compare
            ignore_strand (bool, default False): Whether coordinates must
            be on the same strand to overlap.
        Returns:
            (bool): True if yes, False if no
        """
        strand_match = (
            ignore_strand or self.is_reverse_strand == target.is_reverse_strand
        )
        return (
            self.chrom == target.chrom
            and strand_match
            and not (self.stop < target.start or target.stop < self.start)
        )

    def difference(self, coord, ignore_strand=False):
        """
        Return a list of Onedinate representing the result of removing all
        overlaps of the given Onedinate from this Onedinate.
        Args:
            coord (Onedinate): other Onedinate whose overlaps with this
            Onedinate will be removed from this Onedinate.
            ignore_strand (bool, default False): whether coord must be
            on the same strand as self to return a difference. If True,
            the difference will be returned on the same strand as self
            (if it exists).
        Returns:
            list of Onedinate: a list of Onedinate objects representing
            this Onedinate with overlapping sections of the given Onedinate
            removed.
        """
        if not isinstance(coord, Onedinate):
            raise Exception(
                f"Cannot find difference between self of type {type(self)} and "
                f"coord of type {type(coord)}"
            )
        intersection = self.intersection(coord, ignore_strand=ignore_strand)
        if not intersection:
            return [self]
        else:
            segments = []
            if intersection.start > self.start:
                segments.append(
                    Onedinate(
                        self.chrom,
                        self.start,
                        intersection.start - 1,
                        self.is_reverse_strand,
                    )
                )
            if intersection.stop < self.stop:
                segments.append(
                    Onedinate(
                        self.chrom,
                        intersection.stop + 1,
                        self.stop,
                        self.is_reverse_strand,
                    )
                )
            return segments

    def is_5prime_of(self, coord):
        """
        Return whether the 3' end of self is 5' of the 5' end of coord on the same chromosome/strand.
        """
        if not isinstance(coord, Onedinate):
            raise Exception(
                f"Cannot compare coordinate of type {type(coord)} with self of type {type(self)}"
            )
        elif (
            self.chrom != coord.chrom
            or self.is_reverse_strand != coord.is_reverse_strand
        ):
            return False
        elif self.is_reverse_strand:
            return self.boundary_position_3prime > coord.boundary_position_5prime
        else:
            return self.boundary_position_3prime < coord.boundary_position_5prime

    def is_3prime_of(self, coord):
        """
        Return whether the 5' end of self is 3' of the 3' end of coord on the same chromosome/strand.
        """
        if not isinstance(coord, Onedinate):
            raise Exception(
                f"Cannot compare coordinate of type {type(coord)} with self of type {type(self)}"
            )
        elif (
            self.chrom != coord.chrom
            or self.is_reverse_strand != coord.is_reverse_strand
        ):
            return False
        elif self.is_reverse_strand:
            return self.boundary_position_5prime < coord.boundary_position_3prime
        else:
            return self.boundary_position_5prime > coord.boundary_position_3prime

    def to_zerodinate(self):
        """
        Convert to Zerodinate.
        ATTENSTION: Onedinate cannot indicate positions between two bases like Zerodinate.
        Therefore, chr1:1-1,+ (which indicate the first base of chr1 on + strand) will be converted to chr1:0_1,+
        Returns:
            Zerodinate
        """
        return Zerodinate(self.chrom, self.start - 1, self.stop, self.is_reverse_strand)

    @property
    def length(self):
        return self.stop - self.start + 1

    @property
    def string_value(self):
        """
        Returns: the string representation of a Onedinate, which will be:
        chr1:123-456,+
        chr1:456-789,-
        """
        if self.is_reverse_strand:
            return "{}:{}-{},-".format(self.chrom, self.start, self.stop)
        else:
            return "{}:{}-{},+".format(self.chrom, self.start, self.stop)

    @property
    def string_value_unstranded(self):
        """
        Returns: the string representation of a Onedinate without the strand, which will be:
        chr1:123-456
        chr1:456-789
        """
        return "{}:{}-{}".format(self.chrom, self.start, self.stop)


class Zerodinate(Coordinate):
    def __init__(self, chrom, start=None, stop=None, is_reverse_strand=None):
        super(Zerodinate, self).__init__(chrom, start, stop, is_reverse_strand)

    def is_overlapping_with(self, target, ignore_strand=False):
        """
        Check if the current Zerodinate is overlapping with the target one.
        The strand does affect the result. If ignore_strand is True,
        coordinates must be on the same strand to overlap.
        Args:
            target (Zerodinate): Another Zerodinate to compare
            ignore_strand (bool, default False): Whether coordinates must
            be on the same strand to overlap.
        Returns:
            (bool): True if yes, False if no
        """
        strand_match = (
            ignore_strand or self.is_reverse_strand == target.is_reverse_strand
        )
        return (
            self.chrom == target.chrom
            and strand_match
            and not (self.stop <= target.start or target.stop <= self.start)
        )

    def difference(self, coord, ignore_strand=False):
        """
        Return a list of Zerodinate representing the result of removing all
        overlaps of the given Zerodinate from this Zerodinate.
        Args:
            coord (Zerodinate): other Zerodinate whose overlaps with this
            Zerodinate will be removed from this Zerodinate.
            ignore_strand (bool, default False): whether coord must be
            on the same strand as self to return a difference. If True,
            the difference will be returned on the same strand as self
            (if it exists).
        Returns:
            list of Zerodinate: a list of Zerodinate objects representing
            this Zerodinate with overlapping sections of the given Zerodinate
            removed.
        """
        if not isinstance(coord, Zerodinate):
            raise Exception(
                f"Cannot find difference between self of type {type(self)} and "
                f"coord of type {type(coord)}"
            )
        intersection = self.intersection(coord, ignore_strand=ignore_strand)
        if not intersection:
            return [self]
        else:
            segments = []
            if intersection.start > self.start:
                segments.append(
                    Zerodinate(
                        self.chrom,
                        self.start,
                        intersection.start,
                        self.is_reverse_strand,
                    )
                )
            if intersection.stop < self.stop:
                segments.append(
                    Zerodinate(
                        self.chrom, intersection.stop, self.stop, self.is_reverse_strand
                    )
                )
            return segments

    def is_5prime_of(self, coord):
        """
        Return whether the 3' end of self is 5' of the 5' end of coord on the same chromosome/strand.
        """
        if not isinstance(coord, Zerodinate):
            raise Exception(
                f"Cannot compare coordinate of type {type(coord)} with self of type {type(self)}"
            )
        elif (
            self.chrom != coord.chrom
            or self.is_reverse_strand != coord.is_reverse_strand
        ):
            return False
        elif self.is_reverse_strand:
            return self.boundary_position_3prime >= coord.boundary_position_5prime
        else:
            return self.boundary_position_3prime <= coord.boundary_position_5prime

    def is_3prime_of(self, coord):
        """
        Return whether the 5' end of self is 3' of the 3' end of coord on the same chromosome/strand.
        """
        if not isinstance(coord, Zerodinate):
            raise Exception(
                f"Cannot compare coordinate of type {type(coord)} with self of type {type(self)}"
            )
        elif (
            self.chrom != coord.chrom
            or self.is_reverse_strand != coord.is_reverse_strand
        ):
            return False
        elif self.is_reverse_strand:
            return self.boundary_position_5prime <= coord.boundary_position_3prime
        else:
            return self.boundary_position_5prime >= coord.boundary_position_3prime

    def to_onedinate(self):
        """
        Convert to Onedinate.
        ATTENTION: Onedinate cannot indicate positions between two bases like Zerodinate.
        Therefor, if try to convert Zerodinate chr1:1_1,+ will raise an exception.
        Returns:
            Onedinate
        """
        if self.start == self.stop:
            raise Exception(
                f"Cannot convert {self} to Onedinate. Onedinate cannot represent positions between two bases."
            )
        return Onedinate(self.chrom, self.start + 1, self.stop, self.is_reverse_strand)

    @property
    def length(self):
        return self.stop - self.start

    @property
    def string_value(self):
        """
        Returns: the string representation of a Zerodinate, which will be:
        chr1:123_456,+
        chr1:456_789,-
        """
        if self.is_reverse_strand:
            return "{}:{}_{},-".format(self.chrom, self.start, self.stop)
        else:
            return "{}:{}_{},+".format(self.chrom, self.start, self.stop)

    @property
    def string_value_unstranded(self):
        """
        Returns: the string representation of a Zerodinate without the strand, which will be:
        chr1:123_456
        chr1:456_789
        """
        return "{}:{}_{}".format(self.chrom, self.start, self.stop)
