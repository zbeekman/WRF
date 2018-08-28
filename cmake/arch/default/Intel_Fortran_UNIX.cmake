# Copyright 2018 M. Riechert and D. Meyer. Licensed under the MIT License.

set(checked "-check bounds")
set(preprocess "-fpp")
set(io "-convert big_endian -assume byterecl")
set(optimized "${optimized} -align all -fno-alias")
set(other "-fp-model precise -DNONSTANDARD_SYSTEM_FUNC")