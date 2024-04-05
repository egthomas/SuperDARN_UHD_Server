/***********************************************************************
 * recv_and_hold() function
 * A function to be used to receive samples from the USRP and
 * hold them in the network socket buffer. *rx_data_buffers points to
 * the memory locations of each antenna's samples.   Meant to operate in its
 * own thread context so that it does not block the execution in main()
 * Also handles the former duties of the timing card...
 **********************************************************************/
#include <complex>
#include <vector>
#include <queue>
#include <iostream>
#include <iomanip>
#include <thread>
#include <math.h>

#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/utils/thread_priority.hpp>
#include <uhd/utils/safe_main.hpp>
#include <uhd/utils/static.hpp>
#include <uhd/exception.hpp>
#include <boost/thread.hpp>

#include "rx_worker.h"
#include "usrp_utils.h"
#include "dio.h"

#define DEBUG 1
#ifdef DEBUG
#define DEBUG_PRINT(...) do{ fprintf( stderr, __VA_ARGS__ ); } while( false )
#else
#define DEBUG_PRINT(...) do{ } while ( false )
#endif

#define RX_OFFSET 0


void usrp_rx_worker(
    uhd::usrp::multi_usrp::sptr usrp,
    uhd::rx_streamer::sptr rx_stream,
    std::vector<std::vector<std::complex<int16_t>>> *rx_data_buffer,
    size_t num_requested_samps,
    uhd::time_spec_t start_time,
    int32_t *return_status
){

    DEBUG_PRINT("entering RX_WORKER\n");
 //   fprintf( stderr, "RX WORKER nSamples requested: %i\n", num_requested_samps );
    int nSides = (*rx_data_buffer).size();
    const size_t num_max_request_samps = rx_stream->get_max_num_samps();

 //   fprintf( stderr, "RX WORKER nSides : %i\n", nSides );

    // max number of samples per stream command is limited by UHD
    // to 0x0fffffff samples at the full sampling rate of the ADC (200 MSPS)
    // we request samples after downconversion, so the maximum number of samples that we can request is 0x0fffffff / (200e6 / output_rate)
    // to make things a bit easier, just hardcode a conservative max number of samples
    // also make math easier by chunking at a multiple of the max number of samples per packet..
    //
    
    // then, issue multiple stream commands to gather the full amount of samples

    //setup streaming
    uhd::rx_metadata_t md;
    md.error_code = uhd::rx_metadata_t::ERROR_CODE_NONE;

    double timeout = 5.0;
    
    uhd::stream_cmd_t stream_cmd = uhd::stream_cmd_t::STREAM_MODE_NUM_SAMPS_AND_MORE;
    stream_cmd.num_samps = num_requested_samps/4;
    stream_cmd.stream_now = false;
    stream_cmd.time_spec = offset_time_spec(start_time, RX_OFFSET);

    DEBUG_PRINT("rx_worker requesting %zu samples\n", stream_cmd.num_samps);   
    DEBUG_PRINT("rx_worker requesting %zh (hex) samples\n", stream_cmd.num_samps);   
    DEBUG_PRINT("rx_worker start issue stream\n");   
    usrp->issue_stream_cmd(stream_cmd);

    stream_cmd.stream_now = true;
    usrp->issue_stream_cmd(stream_cmd);
    usrp->issue_stream_cmd(stream_cmd);

    stream_cmd = uhd::stream_cmd_t::STREAM_MODE_NUM_SAMPS_AND_DONE;
    stream_cmd.num_samps = num_requested_samps/4;
    stream_cmd.stream_now = true;

    usrp->issue_stream_cmd(stream_cmd);

    DEBUG_PRINT("rx_worker end issue strem\n");
    
    size_t num_acc_samps = 0;
    std::vector<std::complex<int16_t>*> buff_ptrs(nSides);
/*   for (int iSide=0;iSide<nSides;iSide++) {
      DEBUG_PRINT("side %d \n  ", iSide);
    
      DEBUG_PRINT("  v1: \n  ");
      for (int iSample=0; iSample<10; iSample++)
          DEBUG_PRINT("%i, ", (*rx_data_buffer)[iSide][iSample]);
    }
 */

    DEBUG_PRINT("starting rx_worker while loop\n");
    while(num_acc_samps < num_requested_samps) {
        size_t samp_request = std::min(num_max_request_samps, num_requested_samps - num_acc_samps);
        for (int iSide = 0; iSide < nSides; iSide++) {
            buff_ptrs[iSide] = &((*rx_data_buffer)[iSide][num_acc_samps]);
          //  if (num_acc_samps == 0)
          //     DEBUG_PRINT("rx_worker addr: %p iSide: %d \n    ", buff_ptrs[iSide], iSide);
        }

        size_t num_rx_samps = rx_stream->recv(buff_ptrs , samp_request, md, timeout);
        
        timeout = 0.1;

        //handle the error codes
        if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_TIMEOUT) break;
        if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_LATE_COMMAND) break;
        if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_OVERFLOW) {
           // same error code used for overflow and out of sequence
           if (md.out_of_sequence) {
               std::cerr << "ERROR_CODE_OVERFLOW: out of sequece." << std::endl;
           } else {
               std::cerr << "ERROR_CODE_OVERFLOW: overflow (not out of sequece)" << std::endl;
           }
           break;
        }

        if (md.error_code != uhd::rx_metadata_t::ERROR_CODE_NONE){
            throw std::runtime_error(str(boost::format(
                "Receiver error %s"
            ) % md.strerror()));
        }
                
        //if (DEBUG) {
        //    std::cout << boost::format("Received packet: %u samples, %u full secs, %f frac secs") % num_rx_samps % md.time_spec.get_full_secs() % md.time_spec.get_frac_secs() << std::endl;
        //}
        num_acc_samps += num_rx_samps;
    }

    DEBUG_PRINT("RX_WORKER fetched samples!\n");
    if(DEBUG) std::cout << boost::format("RX_WORKER : %u full secs, %f frac secs") % md.time_spec.get_full_secs() % md.time_spec.get_frac_secs() << std::endl;

    int mute_received_samples = 0;

    if (num_acc_samps != num_requested_samps){
        *return_status=-1;
        uhd::time_spec_t rx_error_time = usrp->get_time_now();
        std::cerr << "Error in receiving samples..(" << rx_error_time.get_real_secs() << ")\t";

        std::cerr << "Error code: " << md.error_code << "\t";
        std::cerr << "Samples rx'ed: " << num_acc_samps << 
            " (expected " << num_requested_samps << ")" << std::endl;
    }

    if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_LATE_COMMAND) {
        uhd::time_spec_t rx_error_time = usrp->get_time_now();
        std::cerr << "Timeout encountered at " << rx_error_time.get_real_secs() << std::endl;
        *return_status=-1;
    }
    if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_TIMEOUT) {
        uhd::time_spec_t rx_error_time = usrp->get_time_now();
        std::cerr << "Timeout encountered at " << rx_error_time.get_real_secs() << std::endl;
        *return_status=-1;
    }
    if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_OVERFLOW){
        uhd::time_spec_t rx_error_time = usrp->get_time_now();
        std::cerr << "Overflow encountered at " << rx_error_time.get_real_secs() << std::endl;

  //      size_t num_rx_samps = rx_stream->recv(buff_ptrs , num_max_request_samps, md, timeout);
  //    std::cerr << "Overflow cleanup: received  " << num_rx_samps << "  , end _of _burst: " << md.end_of_burst <<  std::endl;

/*        while (num_rx_samps != 0 && !md.end_of_burst ) {
            num_rx_samps = rx_stream->recv(buff_ptrs , num_max_request_samps, md, timeout);
            std::cerr << "Overflow cleanup: received  " << num_rx_samps << "  , end _of _burst: " << md.end_of_burst <<  std::endl;
        }
*/
        std::cerr << "Overflow cleanup: finished."   <<  std::endl;

        mute_received_samples = 1;
        *return_status=-1;
    }
    if (md.error_code != uhd::rx_metadata_t::ERROR_CODE_NONE){
        uhd::time_spec_t rx_error_time = usrp->get_time_now();
        std::cerr << "start time: " << start_time.get_real_secs() << std::endl;
        std::cerr << "Unexpected error code " << md.error_code << " encountered at " << rx_error_time.get_real_secs() << std::endl;
        *return_status=-1;
    }
    if (md.out_of_sequence){
        uhd::time_spec_t rx_error_time = usrp->get_time_now();
        std::cerr << "start time: " << start_time.get_real_secs() << std::endl;
        std::cerr << "Packets out of order " << " encountered at " << rx_error_time.get_real_secs() << std::endl;
        *return_status=-1;
    }

    if (mute_received_samples) {
        std::cerr << "Muting received samples because error occurred. " << std::endl;
        for (int iSide = 0; iSide < nSides; iSide++) { // TODO maybe smarter to mute sample while transfer to SHM (and just use special return_status)
           for (int iSample=0; iSample<num_requested_samps; iSample++) {
              (*rx_data_buffer)[iSide][iSample] = 0;
           }
        }
    }


    DEBUG_PRINT("RX_WORKER finished\n");
    return;
}

