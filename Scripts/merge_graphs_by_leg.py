import argparse
from collections import defaultdict
import logging
from pprint import pprint
import sys

from schrodinger.application.desmond import constants
from schrodinger.application.scisol.fep import graph
from schrodinger.utils import log

# local, depends on SCHRODINGER_PYTHONPATH
from pptc.util import mutation as pm

# Configure logging
DEFAULT_LOGGING_LEVEL = logging.WARNING
logger = log.get_output_logger(__name__)
logger.level = DEFAULT_LOGGING_LEVEL


def get_matching_source_edge(dest_edge, title_map, mutation_map):
    '''
    Find and return an edge in `source_graph` that matches `dest_edge`.

    If no matching edge is found, return None.

    Parameters:
    - source_graph (graph.Graph)
    - dest_edge (graph.Edge)
    - title_map (Dict[str, graph.Edge]) mapping `Edge.short_id_title` to Edge
    - mutation_map (Dict[str, graph.Edge]) mapping a mutation string to Edge

    Mutation strings use the same format as protein mutation FEP Node titles,
    e.g. "A-ALA123ASER" (`<chain>-<start_res><resnum>[<inscode>]<end_res>`).

    '''
    logger.debug(f'    dest edge: {dest_edge.short_id_title}')
    dest_pme = pm.ProteinMutationEdge(edge=dest_edge)
    logger.debug(f'    dest mutation: {dest_pme.mutation}')

    mut_match = str(dest_pme.mutation) in mutation_map.keys()
    # mut_is_charged = dest_pme.mutation.is_charged()
    # rev_mut_match = str(reversed(dest_pme.mutation)) in mutation_map.keys()

    source_edge = None
    # Exact title match
    if dest_edge.short_id_title in title_map.keys():
        source_edge = title_map[dest_edge.short_id_title]
    # Mutation match
    elif mut_match:
        source_edge = mutation_map[str(dest_pme.mutation)]
    # Reverse mutation match for neutral mutations
    # elif rev_mut_match and not mut_is_charged:
    #     source_edge = mutation_map[str(reversed(dest_pme.mutation))]
    else:
        logger.warning(
            f"*** WARNING: No match for "
            f"{dest_pme.mutation} "
            f"(edge {dest_edge.short_id_title}) ***"
        )
    return source_edge


def copy_leg(source_edge, dest_edge, source_leg_name, dest_leg_name,
             dg_only=False, wt_min=False):
    '''
    Copy leg information from `source_leg_name` on `source_edge` into
    `dest_leg_name` on `dest_edge`.
    '''
    leg_name_map = {
        'complex': {
            'leg_type': constants.FepLegTypes.COMPLEX,
            'sid_key': graph._KEY_COMPLEX_SID,
        },
        'solvent': {
            'leg_type': constants.FepLegTypes.SOLVENT,
            'sid_key': graph._KEY_SOLVENT_SID,
        },
        'fragment': {
            'leg_type': constants.FepLegTypes.FRAGMENT,
            'sid_key': None,
        },
    }
    dest_leg_type = leg_name_map[dest_leg_name]['leg_type']
    source_leg_dg = source_edge.get_leg_dg_by_name(source_leg_name)
    dest_leg_dg = dest_edge.get_leg_dg_by_name(dest_leg_name)

    sid_key = leg_name_map[source_leg_name]['sid_key']
    if dg_only:
        source_leg_sid = None
    else:
        source_leg_sid = source_edge.get_data(sid_key)

    if wt_min:
        # If this is a direct edge from WT, use the minimum dG value
        n0, _ = source_edge
        is_direct_edge = n0.struc.title == "WT"
        has_dg = source_leg_dg is not None
        print(f'is direct edge: {is_direct_edge} / has_dg: {has_dg}')
        if is_direct_edge and has_dg and (dest_leg_dg.val < source_leg_dg.val):
            logger.warning(
                f"→ NOT copying (dg {source_leg_dg}) from "
                f"{source_edge.short_id_title} ({source_leg_name}) "
                f"to {dest_edge.short_id_title} ({dest_leg_name}) "
                f"because dest dg {dest_leg_dg} is lower"
            )
            return

    # Set the new leg info
    if dest_leg_name in dest_edge.get_leg_names():
        dest_leg = dest_edge.get_leg_by_name(dest_leg_name)
        dest_leg.dg = source_leg_dg
        if sid_key is not None:
            dest_edge.set_data(sid_key, source_leg_sid)
    else:
        dest_edge.add_leg(dest_leg_name, dest_leg_type, dg=source_leg_dg,
                          sid=source_leg_sid)
    logger.info(
        f"→ Copied "
        f"(dg {source_leg_dg}{' + sid' if not dg_only else ''}) "
        f"from {source_edge.short_id_title} ({source_leg_name}) "
        f"to {dest_edge.short_id_title} ({dest_leg_name})"
    )


def delete_leg_from_edge(edge, leg_name):
    '''
    Delete the specified leg info from the given edge.
    '''
    try:
        edge._legs.pop(leg_name)
        logger.debug(f"  - deleted {leg_name} leg from edge {edge.short_id_title}")
    except KeyError:
        logger.debug(f"  - no {leg_name} leg found for edge "
                     f"{edge.short_id_title}, nothing deleted")
        pass


def get_equiv_mutstrs(mutstr, equiv_chains, sep=","):
    '''
    Return a list of possible equivalent mutation strings for the given
    equivalent chains.

    mutstr: mutation string, e.g. "A-ALA123ASER"
    equiv_chains: dict (str -> list of str)
        e.g.  {"AE": ["BC", "CD", "DF"]} means if A is replaced by B, then E
        should be replaced by C, and so forth.  Essentially these are the
        structurally interchangeable complexes.
    '''
    equiv_chain_replacements = [
        list(zip(source_chains, equiv_chains))
        for source_chains, equiv_chains_list in equiv_chains.items()
        for equiv_chains in equiv_chains_list
    ]
    equiv_mutstrs = []
    for chain_replacement_list in equiv_chain_replacements:
        logger.debug(f"chain_replacement_list = {chain_replacement_list}")
        single_mutstr_list = mutstr.split(sep)
        equiv_mutstr_list = [None] * len(single_mutstr_list)
        for i, single_mutstr in enumerate(single_mutstr_list):
            for chain_pair in chain_replacement_list:
                chain_hyphen_pair = [chain + '-' for chain in chain_pair]
                logger.debug(f"checking {single_mutstr} using pair {chain_pair}")
                equiv_mutstr = single_mutstr.replace(*chain_hyphen_pair)
                if equiv_mutstr != single_mutstr:
                    equiv_mutstr_list[i] = equiv_mutstr
                    # only do one replacement per single mutstr per chain pair
                    break
        if None not in equiv_mutstr_list:
            equiv_mutstrs.append(sep.join(equiv_mutstr_list))
    return equiv_mutstrs


def copy_graph_legs(source_graph, dest_graph, source_leg_name, dest_leg_name,
                    dg_only=False, equiv_chains=None, delete_dest_legs=False,
                    wt_min=False):
    '''
    Copy `source_leg_name` from `source_graph` into `dest_leg_name` in `dest_graph`
    for all edges matched between the two graphs.

    If wt_min is True, then for direct Protein FEP edges from the WT node,
    report the minimum leg dG value in the output graph; this reflects an
    assumption that the WT input model easily samples its lowest energy
    conformation in all legs.

    Note: the leg indicated by `source_leg_name` is deleted for all edges in
    `dest_graph`, regardless of whether a matching edge is found in
    `source_graph`.

    '''
    title_to_source_edge_map = {}
    mutation_to_source_edge_map = {}
    for source_edge in source_graph.edges():
        source_pme = pm.ProteinMutationEdge(edge=source_edge)
        logger.debug(f"new source edge: {source_edge.short_id_title}")
        logger.debug(f"  - actual perturbation mutation: {source_pme.mutation}")
        # logger.debug(f"  - bg_mutations: {source_pme.bg_mutations}")

        # Edge titles are always unique.
        title_to_source_edge_map[source_edge.short_id_title] = source_edge

        # FIXME: Edge mutations may not be unique!  Currently we assume all
        # single mutations, with no background mutations in the source graph
        # (i.e. all mutations are unique in the graph, no mutations repeated
        # in different contexts).  However, even if this assumption is not
        # valid for a given graph, it is likely a decent approximation.
        source_mutstr = str(source_pme.mutation)
        mutation_to_source_edge_map[source_mutstr] = source_edge

        # Allow dGs from equivalent chains
        if equiv_chains is not None:
            equiv_mutstrs = get_equiv_mutstrs(source_mutstr, equiv_chains)
            for equiv_mutstr in equiv_mutstrs:
                mutation_to_source_edge_map[equiv_mutstr] = source_edge


    for dest_edge in dest_graph.edges():
        logger.debug(f"\ndest edge: {dest_edge.short_id_title}")
        logger.debug(f"  - legs: {dest_edge.get_leg_names()}")

        if delete_dest_legs:
            # Prevent mixing of `dest_leg_name` data from different graphs
            delete_leg_from_edge(dest_edge, dest_leg_name)

        source_edge = get_matching_source_edge(dest_edge,
                                               title_to_source_edge_map,
                                               mutation_to_source_edge_map)
        if source_edge:
            logger.debug(f"source edge: {source_edge.short_id_title}")
            copy_leg(source_edge, dest_edge, source_leg_name, dest_leg_name,
                     dg_only=dg_only, wt_min=wt_min)


def main(main_fmp, leg_specs, out_fmp=None, fep_type=None, dg_only=False,
         equiv_chains=None, wt_min=False):
    '''
    Main workflow: Copy legs indicated by leg specifications in `leg_specs`
    into a copy of the `main_fmp` graph, and write the resulting output graph.

    '''
    equiv_chains_dict = defaultdict(list)
    if equiv_chains is not None:
        for mapping in equiv_chains:
            source, equiv = mapping.split(':')
            assert len(source) == len(equiv)  # sanity check
            equiv_chains_dict[source].append(equiv)

    # main_graph = graph.Graph.deserialize(main_fmp)
    # dest_graph = copy.deepcopy(main_graph)  # TODO: is this necessary?
    dest_graph = graph.Graph.deserialize(main_fmp)
    for i, leg_spec in enumerate(leg_specs):
        source_fmp, source_leg_name, dest_leg_name = leg_spec.split(":")
        source_graph = graph.Graph.deserialize(source_fmp)
        delete = (i == 1)  # only delete the first time through the dest graph
        copy_graph_legs(source_graph, dest_graph, source_leg_name, dest_leg_name,
                        dg_only=dg_only, equiv_chains=equiv_chains_dict,
                        delete_dest_legs=delete, wt_min=wt_min)

    if fep_type:
        # TODO:  check `fep_type` against supported values
        dest_graph.set_data('fep_type', fep_type)

    logger.info("\nCalculating cycle closure...")
    dest_graph.calc_cycle_closure()

    if out_fmp is None:
        out_fmp = 'merged_by_leg.fmp'
    logger.info("\nWriting output graph...")
    dest_graph.write(out_fmp)
    logger.warning(f"\nWrote {out_fmp}.")
    logger.info("\nDone.")


def parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('fmp',
                        help='primary .fmp graph')
    parser.add_argument('-l', dest='leg_specs', action='append',
                        help='additional graph + leg specification to be '
                             'copied to the primary graph, in the form '
                             '`<fmp_file>:<source_leg_name>:<dest_leg_name>`')
    parser.add_argument('-o', dest='out_fmp', default=None,
                        help='output .fmp file name')
    parser.add_argument('-equiv-chains', dest='equiv_chains', action='append',
                        help='colon-separated mapping of equivalent chains, '
                             'e.g. "AB:CD" for A=C and B=D equivalence')
    parser.add_argument('-fep-type',
                        help='fep_type for output graph')
    parser.add_argument('-dg-only', action='store_true',
                        help='exclude SID report information from output graph')
    parser.add_argument('-wt-min', action='store_true',
                        help='For direct Protein FEP edges from the WT node, '
                        'report the minimum leg dG value in the output graph; '
                        'this reflects an assumption that the WT input model '
                        'easily samples its lowest energy conformation in all '
                        'legs')
    parser.add_argument('-v', dest='verbose', action='store_true',
                        help='increase verbosity')
    parser.add_argument('-debug', action='store_true',
                        help='print debugging info')
    return parser.parse_args(argv)


if __name__ == '__main__':
    args = parse_args(sys.argv[1:])

    # Adjust logging
    if args.debug:
        logger.setLevel(logging.DEBUG)
    elif args.verbose:
        logger.setLevel(logging.INFO)

    main(main_fmp=args.fmp,
         leg_specs=args.leg_specs,
         out_fmp=args.out_fmp,
         fep_type=args.fep_type,
         dg_only=args.dg_only,
         equiv_chains=args.equiv_chains,
         wt_min=args.wt_min)
