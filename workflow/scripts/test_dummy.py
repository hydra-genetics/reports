# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Niklas Mähler"
__copyright__ = "Copyright 2023, Niklas Mähler"
__email__ = "niklas.mahler@regionvasterbotten.se"
__license__ = "GPL-3"


def test_dummy():
    from dummy import dummy
    assert dummy() == 1
