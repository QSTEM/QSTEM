The design of the potential structure is the following:

A general interface that must be adhered to is defined in pot_interface.hpp.  It defines nothing about the data that a potential represents, but only the code interface that outside code uses to manipulate and get data from potential objects.

A general starter class is provided in pot_base.hpp.  This provides the data arrays that are used in the 4 fundamental potential types.  It excludes specific details to how potential gets calculated.

Four subclasses of the class in pot_base.cpp make up the actual concrete classes.  These are 2D/3D potentials, computed either in real-space or using the FFT method.




If you're a developer looking to do something new, first consider subclassing one of the four concrete classes and overriding your desired functionality.  Failing that, go to the CPotential class.  If that's still more than you want, subclass IPotential (you'll be implementing the interface, ensuring the rest of QSTEM knows how to talk to your potential.)

For any new class you want to add, you need to add it to the factory to make it generally available.  In your new class' header file, paste a section like this:

private:
	friend class CPotFactory;
	// Create an instance of this class, wrapped in a shared ptr
	//     This should not be inherited - any subclass needs its own implementation.
	static PotPtr __stdcall Create(){return PotPtr(new C2DPotential());}

(replace C2DPotential with the name of your class)

In pot_factory.cpp, add a new include line for your file.  In the CPotFactory constructor, call the Register method with a string identifier for your class, followed by a pointer to the Create method you called.  Just copy/pase one of the existing lines, and change the class name to your class name.