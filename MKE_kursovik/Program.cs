using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;

namespace MKE_kursovik
{

    public class Program
    {
        IO io;
        Solve solve;

        static void Main(string[] args)
        {
            Program program = new Program();

			//program.io = new IO(5);
			//program.solve = new Solve(program.io, 5);

			//program.io = new IO(1);
			//program.solve = new(program.io, 1);

			program.io = new IO(2);
			program.solve = new(program.io, 2);

			//program.io = new IO(3);
			//program.solve = new(program.io, 3);

			//program.io = new IO(4);
			//program.solve = new(program.io, 4);

		}
    }
}
