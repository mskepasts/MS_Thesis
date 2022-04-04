using System;
using System.Drawing;
using Grasshopper;
using Grasshopper.Kernel;

namespace oneNode
{
    public class oneNodeInfo : GH_AssemblyInfo
    {
        public override string Name => "oneNode";

        //Return a 24x24 pixel bitmap to represent this GHA library.
        public override Bitmap Icon => null;

        //Return a short string describing the purpose of this GHA library.
        public override string Description => "";

        public override Guid Id => new Guid("9EC875D9-2571-458E-923A-7851BB33E5D4");

        //Return a string identifying you or your company.
        public override string AuthorName => "";

        //Return a string representing your preferred contact details.
        public override string AuthorContact => "";
    }
}