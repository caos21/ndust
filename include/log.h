// Based on
/*
 *          Copyright Andrey Semashev 2007 - 2015.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef LOG_H
#define LOG_H

#include <iostream>
// #include <memory>

#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/support/date_time.hpp>
#include <boost/log/sinks/text_file_backend.hpp>
#include <boost/thread/thread.hpp>

namespace logging = boost::log;
namespace src = boost::log::sources;
namespace expr = boost::log::expressions;
namespace keywords = boost::log::keywords;
namespace sinks = boost::log::sinks;
using namespace logging::trivial;

// typedef sinks::synchronous_sink< sinks::text_file_backend > sink_t;

namespace blog {
  
//   boost::shared_ptr< sink_t > g_file_sink;
//   boost::shared_ptr< sinks::synchronous_sink< sinks::text_file_backend > > pLogSink;
//   shared_ptr< synchronous_sink< text_file_backend > > pLogSink;
  
//   inline
//   void init(std::string logname) {
//     boost::shared_ptr< sinks::synchronous_sink< sinks::text_file_backend > > pLogSink;
//     pLogSink =
//     logging::add_file_log (logname+".log");
// //     std::shared_ptr< sinks::synchronous_sink< sinks::text_file_backend > > pLogSink =
// //     logging::add_file_log (
// //       keywords::file_name = logname+".log",
// //       // This makes the sink to write log records that look like this:
// //       // YYYY-MM-DD HH:MI:SS: [normal] A normal severity message, etc
// //       keywords::format = (
// //         expr::stream
// //             << expr::format_date_time< boost::posix_time::ptime >("TimeStamp", "%Y-%m-%d %H:%M:%S")
// //             << ": [" << severity
// //             << "] " << expr::smessage
// //       )
// //     );
//   }
//   inline
//   void stop() {
//     logging::core::get()->remove_sink(g_file_sink);
//     g_file_sink.reset();
//   }
//   static boost::shared_ptr< sinks::synchronous_sink< sinks::text_file_backend > > pLogSink;
// //   shared_ptr< synchronous_sink< text_file_backend > > pLogSink;
//
   inline
  void init(std::string logname) {
    auto sink = logging::add_file_log (
      keywords::file_name = logname+".log",
      keywords::auto_flush = true,//Write after each call
      // This makes the sink to write log records that look like this:
      // YYYY-MM-DD HH:MI:SS: [normal] A normal severity message, etc
      keywords::format = (
        expr::stream
            << expr::format_date_time< boost::posix_time::ptime >("TimeStamp", "%Y-%m-%d %H:%M:%S")
            << ": [" << severity
            << "] " << expr::smessage
      )
    );
  }
// //   inline
// //   void stop() {
// //     logging::core::get()->remove_sink(g_file_sink);
// //     g_file_sink.reset();
// //   }
}

#endif
