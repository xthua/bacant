# -*- coding: utf-8 -*-
import os
from subprocess import call
import colorlog

from . import EmptyFileError

_log = colorlog.getLogger(__name__)


def find_integrase(replicon_id, prot_file, out_dir, cfg):
    """
    Call Prodigal for Gene annotation and hmmer to find integrase, either with phage_int
    HMM profile or with intI profile.

    :param str replicon_id: The Replicon identifier to search integrase into
    :param str prot_file: the path to the fasta file containing the translation of the replicon.
    :param str out_dir: the relative path to the directory where prodigal outputs will be stored
    :param cfg: the configuration
    :type cfg: a :class:`integron_finder.config.Config` object
    :returns: None, the results are written on the disk
    """

    intI_hmm_out = os.path.join(out_dir, replicon_id + "_intI.res")
    hmm_cmd = []
    if not os.path.exists(prot_file):
        msg = "The protein file: '{}' does not exists cannot perform hmmsearch on it.".format(prot_file)
        _log.warning(msg)
        raise RuntimeError(msg)
    elif os.path.getsize(prot_file) == 0:
        msg = "The protein file: '{}' is empty cannot perform hmmsearch on it.".format(prot_file)
        _log.warning(msg)
        raise EmptyFileError(msg)

    if not os.path.isfile(intI_hmm_out):
        hmm_cmd.append([cfg.hmmsearch,
                        "--cpu", str(cfg.cpu),
                        "--tblout", os.path.join(out_dir, replicon_id + "_intI_table.res"),
                        "-o", intI_hmm_out,
                        cfg.model_integrase,
                        prot_file])

    phage_hmm_out = os.path.join(out_dir, replicon_id + "_phage_int.res")
    if not os.path.isfile(phage_hmm_out):
        hmm_cmd.append([cfg.hmmsearch,
                        "--cpu", str(cfg.cpu),
                        "--tblout", os.path.join(out_dir, replicon_id + "_phage_int_table.res"),
                        "-o", phage_hmm_out,
                        cfg.model_phage_int,
                        prot_file])

    for cmd in hmm_cmd:
        try:
            _log.debug("run hmmsearch: {}".format(' '.join(cmd)))
            returncode = call(cmd)
        except Exception as err:
            raise RuntimeError("{0} failed : {1}".format(' '.join(cmd), err))
        if returncode != 0:
            raise RuntimeError("{0} failed return code = {1}".format(' '.join(cmd), returncode))
