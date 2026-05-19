set(B2bLIB_SOURCES
    convgpx.c
    convbin.c
    convkml.c
    convrnx.c
    datum.c
    download.c
    ephemeris.c
    geoid.c
    gis.c
    ionex.c
    lambda.c
    options.c
    pntpos.c
    postpos.c
    ppp.c
    ppp_ar.c
    ppp_ar_passbypass.c
    ppp_ar_integration.c
    preceph.c
    rcvraw.c
    rinex.c
    rtcm2.c
    rtcm3.c
    rtcm3e.c
    rtcm.c
    rtkcmn.c
    rtkpos.c
    rtkrcv.c
    rtksvr.c
    sbas.c
    solution.c
    stream.c
    streamsvr.c
    tides.c
    tle.c
    trace.c
    vt.c
    B2b.c
    # reciver decoder
    rcv/binex.c
    rcv/crescent.c
    rcv/javad.c
    rcv/novatel.c
    rcv/nvs.c
    rcv/rt17.c
    rcv/septentrio.c
    rcv/skytraq.c
    rcv/swiftnav.c
    rcv/tersus.c
    rcv/ublox.c
    rcv/unicore.c
    rcv/sinan.c
    # fortran to c
    f2c/cal2jd.c
    f2c/dat.c
    f2c/dehanttideinel.c
    f2c/gmf.c
    f2c/gpt.c
    f2c/norm8.c
    f2c/sprod.c
    f2c/st1idiu.c
    f2c/st1isem.c
    f2c/st1l1.c
    f2c/step2diu.c
    f2c/step2lon.c
    f2c/vmf1.c
    f2c/vmf1_ht.c
    f2c/zero_vec8.c
)

add_library(B2bLib STATIC ${B2bLIB_SOURCES})

target_link_libraries(B2bLib PUBLIC m)

if(UNIX)
    target_link_libraries(B2bLib PUBLIC pthread)
elseif(WIN32)
    target_compile_options(B2bLib PUBLIC -mthreads)
    target_link_libraries(B2bLib PUBLIC ws2_32 winmm)
    message(STATUS "MinGW library path: ${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}")
endif()

target_include_directories(B2bLib PUBLIC ${PROJECT_SOURCE_DIR}/include)

target_compile_definitions(B2bLib PUBLIC
    ENAGLO
    ENAGAL
    ENACMP
    ENAQZS
    # LAPACK
)
