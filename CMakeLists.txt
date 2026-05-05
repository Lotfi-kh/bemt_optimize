px4_add_library(bemt_model
	bemt_model.cpp
)

px4_add_module(
	MODULE modules__bemt_optimizer
	MAIN bemt_optimize
	SRCS
		bemt_optimizer.cpp
	DEPENDS
		bemt_model
)

px4_add_unit_gtest(SRC bemt_model_test.cpp LINKLIBS bemt_model)
