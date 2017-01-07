// Copyright (c) 2016, the SDSL Project Authors.  All rights reserved.
// Please see the AUTHORS file for details.  Use of this source code is governed
// by a BSD license that can be found in the LICENSE file.
#ifndef INCLUDED_SDSL_RAM_FSTREAMBUF
#define INCLUDED_SDSL_RAM_FSTREAMBUF

#include <fstream>
#include <vector>
#include "ram_fs.hpp"

namespace sdsl {

class ram_filebuf : public std::streambuf
{
    private:
        // TODO:  also store filename/descriptor to implement is_open ???
        ram_fs::content_type* m_ram_file = nullptr;  // file handle
        void pbump64(std::ptrdiff_t);

public:
	virtual ~ram_filebuf();

	ram_filebuf();
	ram_filebuf(ram_fs::content_type& ram_file);

	std::streambuf* open(const std::string s, std::ios_base::openmode mode);

	bool is_open();

	ram_filebuf* close();

	pos_type seekpos(pos_type				 sp,
					 std::ios_base::openmode which = std::ios_base::in |
													 std::ios_base::out) override;

	pos_type pubseekoff(off_type				off,
						std::ios_base::seekdir  way,
						std::ios_base::openmode which = std::ios_base::in | std::ios_base::out);

	pos_type pubseekpos(pos_type				sp,
						std::ios_base::openmode which = std::ios_base::in | std::ios_base::out);


	std::streamsize xsputn(const char_type* s, std::streamsize n) override;

	int sync() override;

	int_type overflow(int_type c = traits_type::eof()) override;
};
}

#endif
