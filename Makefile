all : HAL2 HAL2-P HAS2 HAS2-P

HAL2 : RAL_new.h alignment.h charSet.h matrix.h tool_box.h RAL_new.cpp alignment.cpp charSet.cpp main_BU_new1.cpp matrix.cpp optimization.h optimization.cpp core.h core.cpp rateMatrixResult.h rateMatrixResult.cpp parameters.h parameters.cpp user_options.h user_options.cpp lbfgsb_new.h lbfgsb_new.cpp definitions.h gradient.h gradient.cpp jacobi_eigenvalue.h jacobi_eigenvalue.cpp tool_box.cpp RAL_common.h RAL_common.cpp global.h global.cpp simd.h
	g++ -o HAL2 RAL_new.cpp alignment.cpp charSet.cpp main_BU_new1.cpp matrix.cpp optimization.cpp core.cpp rateMatrixResult.cpp parameters.cpp user_options.cpp lbfgsb_new.cpp gradient.cpp jacobi_eigenvalue.cpp tool_box.cpp RAL_common.cpp global.cpp -lpthread -O3 -msse4.2 

HAS2 : alignment.h charSet.h matrix.h tool_box.h alignment.cpp charSet.cpp main_RAS.cpp matrix.cpp optimization.h optimization.cpp core.h core.cpp rateMatrixResult.h rateMatrixResult.cpp parameters.h parameters.cpp RAS.h RAS.cpp user_options.h user_options.cpp lbfgsb_new.h lbfgsb_new.cpp definitions.h gradient.h gradient.cpp jacobi_eigenvalue.h jacobi_eigenvalue.cpp tool_box.cpp RAL_common.h RAL_common.cpp global.h global.cpp simd.h
	g++ -o HAS2 alignment.cpp charSet.cpp main_RAS.cpp matrix.cpp optimization.cpp core.cpp rateMatrixResult.cpp parameters.cpp RAS.cpp user_options.cpp lbfgsb_new.cpp gradient.cpp jacobi_eigenvalue.cpp tool_box.cpp RAL_common.cpp global.cpp -lpthread -O3 -msse4.2

HAL2-P : RAL_new.h alignment.h charSet.h matrix.h tool_box.h RAL_new.cpp alignment.cpp charSet.cpp main_BU_new1.cpp matrix.cpp optimization.h optimization.cpp core.h core.cpp rateMatrixResult.h rateMatrixResult.cpp parameters.h parameters.cpp user_options.h user_options.cpp lbfgsb_new.h lbfgsb_new.cpp definitions.h gradient.h gradient.cpp jacobi_eigenvalue.h jacobi_eigenvalue.cpp tool_box.cpp RAL_common.h RAL_common.cpp global.h global.cpp simd.h
	g++ -DC_REG=1 -o HAL2-P RAL_new.cpp alignment.cpp charSet.cpp main_BU_new1.cpp matrix.cpp optimization.cpp core.cpp rateMatrixResult.cpp parameters.cpp user_options.cpp lbfgsb_new.cpp gradient.cpp jacobi_eigenvalue.cpp tool_box.cpp RAL_common.cpp global.cpp -lpthread -O3 -msse4.2

HAS2-P : alignment.h charSet.h matrix.h tool_box.h alignment.cpp charSet.cpp main_RAS.cpp matrix.cpp optimization.h optimization.cpp core.h core.cpp rateMatrixResult.h rateMatrixResult.cpp parameters.h parameters.cpp RAS.h RAS.cpp user_options.h user_options.cpp lbfgsb_new.h lbfgsb_new.cpp definitions.h gradient.h gradient.cpp jacobi_eigenvalue.h jacobi_eigenvalue.cpp tool_box.cpp RAL_common.h RAL_common.cpp global.h global.cpp simd.h
	g++ -DC_REG=1 -o HAS2-P alignment.cpp charSet.cpp main_RAS.cpp matrix.cpp optimization.cpp core.cpp rateMatrixResult.cpp parameters.cpp RAS.cpp user_options.cpp lbfgsb_new.cpp gradient.cpp jacobi_eigenvalue.cpp tool_box.cpp RAL_common.cpp global.cpp -lpthread -O3 -msse4.2

clean :
	rm -f HAL2 HAL2-P HAS2 HAS2-P
