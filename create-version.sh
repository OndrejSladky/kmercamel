#!/usr/bin/env bash
cat <(echo "#define VERSION \"") src/version <(echo -) <(git rev-parse HEAD | head -c6) <(echo \") | tr -d '\n' > src/version.h
