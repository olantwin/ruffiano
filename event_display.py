import matplotlib.pyplot as plt
import ROOT
import argparse


DETECTORS = [
"target_1",
"ScoringPlane1_box_1",
"ScoringPlane2_box_1",
"ScoringPlane3_box_1",
"ScoringPlane4_box_1",
]

def main():
    """Run simple event display."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "inputfile",
        help="""Simulation results to use as input. """
        """Supports retrieving file from EOS via the XRootD protocol.""",
    )
    parser.add_argument(
        "-g",
        "--geofile",
        help="""Geometry file to use as input. """
        """Supports retrieving file from EOS via the XRootD protocol.""",
    )
    parser.add_argument(
        "--plots", help="Make nice plots as pdf and png", action="store_true"
    )
    args = parser.parse_args()
    with ROOT.TFile.Open(args.inputfile) as f, ROOT.TFile.Open(args.geofile) as g:
        tree = f["cbmsim"]
        
        for i, event in enumerate(tree):
            position_dict = {}
            for hit in event.sco1_Point:
                track_id = hit.GetTrackID()
                if track_id == -2:
                    continue
                point = ROOT.Math.XYZPoint(hit.GetX(), hit.GetY(), hit.GetZ())
                if track_id in position_dict:
                    position_dict[track_id].append(point)
                else:
                    position_dict[track_id] = [
                        point,
                    ]
            for hit in event.sco2_Point:
                track_id = hit.GetTrackID()
                if track_id == -2:
                    continue
                point = ROOT.Math.XYZPoint(hit.GetX(), hit.GetY(), hit.GetZ())
                if track_id in position_dict:
                    position_dict[track_id].append(point)
                else:
                    position_dict[track_id] = [
                        point,
                    ]
            for hit in event.sco3_Point:
                track_id = hit.GetTrackID()
                if track_id == -2:
                    continue
                point = ROOT.Math.XYZPoint(hit.GetX(), hit.GetY(), hit.GetZ())
                if track_id in position_dict:
                    position_dict[track_id].append(point)
                else:
                    position_dict[track_id] = [
                        point,
                    ]
            for hit in event.sco4_Point:
                track_id = hit.GetTrackID()
                if track_id == -2:
                    continue
                point = ROOT.Math.XYZPoint(hit.GetX(), hit.GetY(), hit.GetZ())
                if track_id in position_dict:
                    position_dict[track_id].append(point)
                else:
                    position_dict[track_id] = [
                        point,
                    ]
            fig, axes = plt.subplots(2, 1, sharex="col")
            plt.tight_layout(h_pad=0, w_pad=0)
            for i, track in enumerate(event.MCTrack):
                start_point = ROOT.Math.XYZPoint(
                    track.GetStartX(), track.GetStartY(), track.GetStartZ()
                )
                if i in position_dict:
                    position_dict[i].insert(0, start_point)
                else:
                    position_dict[i] = [
                        start_point,
                    ]
                positions = [
                    [point.x(), point.y(), point.z()] for point in position_dict[i]
                ]
                positions = list(zip(*positions))
                axes[0].plot(positions[2], positions[0], marker="x")
                axes[1].plot(positions[2], positions[1], marker="x")
            axes[0].set_xlim(-100, 800)
            axes[0].set_ylim(-200, +200)
            axes[1].set_ylim(-200, +200)
            plt.show()
            if args.plots:
                fig.savefig(f"event_display_{i}.png")
                fig.savefig(f"event_display_{i}.pdf")
            break


if __name__ == "__main__":
    ROOT.gErrorIgnoreLevel = ROOT.kWarning
    ROOT.gROOT.SetBatch(True)
    main()
