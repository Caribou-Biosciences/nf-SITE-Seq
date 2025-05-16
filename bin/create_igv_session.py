#!/usr/bin/env python3

############################################################################
############################################################################
## Create an IGV session file (XML format) containing panels/track for    ##
## the input genome FASTA, BAMs, and BigWigs files                        ##
############################################################################
############################################################################


import os
import argparse
from typing import List
import xml.etree.ElementTree as ET


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome-ref", type=str, required=True)
    parser.add_argument(
        "--panel",
        dest="panels",
        required=False,
        nargs="+",
        action="append",
    )
    parser.add_argument("--out-path", type=str, required=True)
    args = parser.parse_args()
    session = IGVSession(args.genome_ref)
    session.write_xml(args.out_path, args.panels)


class IGVSession:
    def __init__(self, genome_path: str) -> None:
        self.genome_path = genome_path

    @staticmethod
    def indent_xml(elem: ET.Element, level: int = 0) -> None:
        i = "\n" + level * "  "
        if len(elem):
            if not elem.text or not elem.text.strip():
                elem.text = i + "  "
            if not elem.tail or not elem.tail.strip():
                elem.tail = i
            for elem in elem:
                IGVSession.indent_xml(elem, level + 1)
            if not elem.tail or not elem.tail.strip():
                elem.tail = i
        else:
            if level and (not elem.tail or not elem.tail.strip()):
                elem.tail = i

    def add_ref_seq_track(self, panel: ET.SubElement) -> None:
        ET.SubElement(
            panel,
            "Track",
            attributeKey="Reference sequence",
            clazz="org.broad.igv.track.SequenceTrack",
            fontSize="10",
            id="Reference sequence",
            sequenceTranslationStrandValue="+",
            shouldShowTranslation="false",
            visible="true",
        )

    def add_bam_track(self, panel: ET.SubElement, bam_path: str) -> None:
        file_name = os.path.basename(bam_path)
        coverage_track = ET.SubElement(
            panel,
            "Track",
            attributeKey=f"{file_name} Coverage",
            autoScale="true",
            clazz="org.broad.igv.sam.CoverageTrack",
            id=f"{bam_path}_coverage",
            name=f"{file_name} Coverage",
            visible="true",
        )
        ET.SubElement(
            coverage_track,
            "DataRange",
            baseline="0.0",
            drawBaseline="true",
            flipAxis="false",
            minimum="0.0",
            maximum="100.0",
            type="LINEAR",
        )
        alignment_track = ET.SubElement(
            panel,
            "Track",
            attributeKey=file_name,
            clazz="org.broad.igv.sam.AlignmentTrack",
            id=bam_path,
            name=file_name,
            color="185,185,185",
            displayMode="COLLAPSED",
            groupByStrand="false",
            visible="true",
        )
        ET.SubElement(
            panel,
            "Track",
            attributeKey=f"{file_name} Junctions",
            autoScale="false",
            clazz="org.broad.igv.sam.SpliceJunctionTrack",
            id=bam_path.replace(".bam", ".bam_junctions"),
            name=f"{file_name} Junctions",
            groupByStrand="false",
            visible="false",
        )
        ET.SubElement(
            alignment_track,
            "RenderOptions",
            colorOption="READ_STRAND",
        )

    def add_bigwig_track(self, panel: ET.SubElement, bw_path: str) -> None:
        file_name = os.path.basename(bw_path)
        bigwig_track = ET.SubElement(
            panel,
            "Track",
            attributeKey=file_name,
            autoScale="true",
            clazz="org.broad.igv.track.DataSourceTrack",
            id=bw_path,
            name=file_name,
            renderer="BAR_CHART",
            windowFunction="none",
            visible="true",
        )
        ET.SubElement(
            bigwig_track,
            "DataRange",
            baseline="0.0",
            drawBaseline="true",
            flipAxis="false",
            minimum="0.0",
            type="LINEAR",
        )

    def write_xml(self, out_path: str, panels: List[List[str]]) -> None:
        session = ET.Element(
            "Session", genome=self.genome_path, locus="All", version="8"
        )
        resources = ET.SubElement(session, "Resources")

        for panel_idx, panel_paths in enumerate(panels):
            panel = ET.SubElement(session, "Panel", name=f"Panel{panel_idx}")
            for data_path in panel_paths:
                if data_path.endswith(".bam"):
                    resource_type = "bam"
                    self.add_bam_track(panel, data_path)
                elif data_path.endswith(".bw"):
                    resource_type = "bw"
                    self.add_bigwig_track(panel, data_path)
                else:
                    raise ValueError(f"Unrecognized file type: {data_path}")

                ET.SubElement(resources, "Resource", path=data_path, type=resource_type)

        feature_panel = ET.SubElement(session, "Panel", name="FeaturePanel")
        self.add_ref_seq_track(feature_panel)

        hidden_attributes = ET.SubElement(session, "HiddenAttributes")
        ET.SubElement(hidden_attributes, "Attribute", name="DATA FILE")
        ET.SubElement(hidden_attributes, "Attribute", name="DATA TYPE")
        ET.SubElement(hidden_attributes, "Attribute", name="NAME")

        IGVSession.indent_xml(session)
        et = ET.ElementTree(session)
        et.write(out_path)


if __name__ == "__main__":
    main()
