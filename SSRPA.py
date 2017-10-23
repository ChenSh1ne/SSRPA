#!/usr/bin/env python
import os
from ssrlib import commands


if __name__ == '__main__':
    args = commands.parse_args()
    args.func(args)



