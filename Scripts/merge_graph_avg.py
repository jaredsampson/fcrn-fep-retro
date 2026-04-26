from argparse import ArgumentParser
from collections import defaultdict
from pprint import pprint, pformat
from typing import List
from typing import Optional
from typing import Union

from schrodinger.application.desmond.constants import PROTEIN_FEP_TYPES
from schrodinger.application.desmond.measurement import Measurement
from schrodinger.application.scisol.fep import fmpdb
from schrodinger.application.scisol.fep import graph
from schrodinger.application.scisol.fep.cmd.reanalyze_fep \
    import _get_fep_legs_from_fep_type
from schrodinger.application.scisol.fep.utils import check_graphs_for_merging


def merge(g0: "graph.Graph", g: "graph.Graph", check_compatibility=True):
    if check_compatibility:
        msg = check_graphs_for_merging(g0, g)
        if msg:
            raise RuntimeError(msg)
    g0.merge(g)


def update_edge_to_leg_dgs_map(
    ggraph: graph.Graph,
    edge_to_leg_dgs_map: Optional[dict] = None
) -> dict:
    '''
    Update the `edge_to_leg_dgs_map` with leg dG data from `ggraph`.
    '''
    fep_legs = _get_fep_legs_from_fep_type(ggraph.fep_type)
    if edge_to_leg_dgs_map is None:
        edge_to_leg_dgs_map = {}
    for edge in ggraph.edges():
        if edge.id not in edge_to_leg_dgs_map:
            edge_to_leg_dgs_map[edge.id] = defaultdict(list)
        for leg in fep_legs:
            dg = edge.get_leg_dg_by_name(leg)
            if dg is not None:
                edge_to_leg_dgs_map[edge.id][leg].append(dg)
    return edge_to_leg_dgs_map


def store_average_dgs(ggraph: graph.Graph, edge_to_leg_dgs_map: dict):
    '''
    Store average dG values for each leg in each edge.
    '''
    print("\nStoring average leg dG values for all edges"
          "\n===========================================")
    for edge in ggraph.edges():
        print(f"Processing edge {edge.short_id_title}:")
        for leg_name in edge_to_leg_dgs_map[edge.id]:
            leg = edge.get_leg_by_name(leg_name)
            if leg is None:
                continue
            if not edge_to_leg_dgs_map[edge.id][leg_name]:
                leg.dg = None
            dgs = edge_to_leg_dgs_map[edge.id][leg_name]
            avg_dg = sum(dgs) / len(dgs)
            leg.dg = avg_dg

            # Logging
            dg_strs = [str(dg) for dg in dgs]
            all_dgs_str = dg_strs[0] if len(dgs) == 1 else pformat(dg_strs, indent=6)
            print(f"  {leg_name}:")
            print("    - leg dGs:", all_dgs_str)
            print("    - average dG:", avg_dg)


def store_min_dgs_for_direct_edges_from_wt(ggraph: graph.Graph,
                                           edge_to_leg_dgs_map: dict):
    '''
    Store leg dG values for direct edges from WT to mutant that minimize the
    energy of the mutant relative to the WT for each leg.
    '''
    print("\nStoring minimum leg dG values for direct edges from WT"
          "\n======================================================")
    for edge in ggraph.edges():
        print(f"Processing edge {edge.short_id_title}:")
        # Check if this is a direct
        t0, t1 = [n.struc.title for n in edge]
        if t0 == "WT":
            func = min
            print("  Edge is direct from WT; using minimum leg dG")
        elif t1 == "WT":
            # In case the direction is somehow reversed, we still want to
            # store the dG value that minimizes the energy of the mutant
            # relative to the WT.
            func = max
            print("  Note: edge direction is reversed; using maximum leg dG")
        else:
            continue

        for leg_name in edge_to_leg_dgs_map[edge.id]:
            leg = edge.get_leg_by_name(leg_name)
            if leg is None:
                continue
            if not edge_to_leg_dgs_map[edge.id][leg_name]:
                leg.dg = None
            else:
                # Get the index of the leg dG value (from the .val attribute of the Measurement object)
                # that minimizes the energy of the mutant relative to the WT using the provided function.
                measurements = edge_to_leg_dgs_map[edge.id][leg_name]
                vals = [m.val for m in measurements]
                uncs = [m.unc for m in measurements]

                # get the index (or indices) of the value identified by func
                keep_val = func(vals)
                matching_indices = [
                    index
                    for index, val in enumerate(vals)
                    if val == keep_val
                ]

                # If there are somehow multiple matching indices (should be
                # very rare), we choose the one with the smallest uncertainty.
                keep_idx = sorted(matching_indices, key=lambda i: uncs[i])[0]
                leg.dg = Measurement(vals[keep_idx], uncs[keep_idx])
                print(f"  {leg_name}:")
                print(f"    - measurements:", measurements)
                print(f"    - kept value at index {keep_idx}:", leg.dg)


def main(argv=None):
    # type: Union(List[str], None)
    usage = """
    $SCHRODINGER/run -FROM scisol merge_graph.py graph1.fmp graph2.fmp -o graph.fmp

    Merge graph in file graph1.fmp and graph2.fmp into one graph,
    and write output as graph.fmp.
    """
    parser = ArgumentParser(
        usage=usage, description='Merge a set of graphs into a single graph.')
    parser.add_argument("fmps", help="input file names", nargs='+')
    parser.add_argument("-sc",
                        "-skip-compatibility-checking",
                        help="If given, skip checking graph "
                        "compatibility",
                        dest="skip_graph_checking",
                        action="store_true")
    parser.add_argument(
        "-fmpdb",
        dest="fmpdb",
        metavar="FMPDB",
        default=None,
        help="Update the path to the associated fmpdb with the given <FMPDB>")
    parser.add_argument("-o",
                        dest="fmp_out",
                        default="merged_out.fmp",
                        help="output file name, defaults to 'merged_out.fmp'")
    parser.add_argument('-average', action='store_true',
                        help='For repeated edges, report average leg dG values '
                        'in the output graph')
    parser.add_argument('-wt-min', action='store_true',
                        help='For direct Protein FEP edges from the WT node, '
                        'report the minimum leg dG value in the output graph; '
                        'this reflects an assumption that the WT input model '
                        'easily samples its lowest energy conformation in all '
                        'legs')
    args = parser.parse_args(argv)

    g0 = graph.Graph.deserialize(args.fmps[0])
    edge_to_leg_dgs_map = update_edge_to_leg_dgs_map(g0)

    for fn in args.fmps[1:]:
        g = graph.Graph.deserialize(fn)
        edge_to_leg_dgs_map = update_edge_to_leg_dgs_map(g, edge_to_leg_dgs_map)
        merge(g0, g, check_compatibility=not args.skip_graph_checking)

    if args.average:
        store_average_dgs(g0, edge_to_leg_dgs_map)
    if args.wt_min:
        if g0.fep_type not in PROTEIN_FEP_TYPES:
            raise RuntimeError(
                "The -wt-min option is only applicable to Protein FEP graphs.")
        # For direct edges from WT, overwrite the average
        store_min_dgs_for_direct_edges_from_wt(g0, edge_to_leg_dgs_map)

    if args.fmpdb is not None:
        g0.fmpdb = fmpdb.FmpDbReader(args.fmpdb)
    g0.calc_cycle_closure()
    g0.write(args.fmp_out)


if "__main__" == __name__:
    main()
