#!/usr/bin/env bash
cat <(echo "#define VERSION \"") src/version <(echo -) <(git rev-parse HEAD) <(echo \") | tr -d '\n' > src/version.h
