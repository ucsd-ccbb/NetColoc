# -*- coding: utf-8 -*-

import click

from .netprop_zscore import *


@click.command(context_settings={"ignore_unknown_options": True})
@click.argument('seed-gene-file', type=click.Path(exists=True, resolve_path=True))
@click.option('--seed-gene-file-delimiter')
@click.option('--num-reps', default=10, type=int)
@click.option('--alpha', '-a', default=0.5, type=float)
@click.option('--minimum-bin-size', default=10, type=int)
@click.option('--interactome-file')
@click.option('--interactome-uuid', default='f93f402c-86d4-11e7-a10d-0ac135e8bacf', type=click.UUID)
@click.option('--ndex-server', default='ndexbio.org')
@click.option('--ndex-user')
@click.option('--ndex-password')
@click.option('--out-name', default='out')
@click.option('--save-z-scores/--no-save-z-scores', default=True)
@click.option('--save-final-heat/--no-save-final-heat', default=False)
@click.option('--save-random-final-heats/--no-save-random-final-heats', default=False)
def main(seed_gene_file, seed_gene_file_delimiter, num_reps, alpha,
         minimum_bin_size, interactome_file,
         interactome_uuid, ndex_server, ndex_user, ndex_password,
         out_name, save_z_scores, save_final_heat, save_random_final_heats):
    """

    :param seed_gene_file:
    :param seed_gene_file_delimiter:
    :param num_reps:
    :param alpha:
    :param minimum_bin_size:
    :param interactome_file:
    :param interactome_uuid:
    :param ndex_server:
    :param ndex_user:
    :param ndex_password:
    :param out_name:
    :param save_z_scores:
    :param save_final_heat:
    :param save_random_final_heats:
    :return:
    """
    netprop_zscore(seed_gene_file, seed_gene_file_delimiter=seed_gene_file_delimiter,
                   num_reps=num_reps, alpha=alpha, minimum_bin_size=minimum_bin_size,
                   interactome_file=interactome_file, interactome_uuid=str(interactome_uuid),
                   ndex_server=ndex_server, ndex_user=ndex_user, ndex_password=ndex_password,
                   out_name=out_name, save_z_scores=save_z_scores,
                   save_final_heat=save_final_heat,
                   save_random_final_heats=save_random_final_heats)


if __name__ == "__main__":  # pragma: no cover
    main()
