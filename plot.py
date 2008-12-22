from ariadne import Box
import Tkinter as Tk

def plot(boxes):
    root = Tk.Tk()
    w = Tk.Canvas(root, width=500, height=500)
    w.pack()

    w.width=500
    w.height=500
    w.left=-1.
    w.right=2.
    w.bottom=-1.
    w.top=2.

    for box in boxes:
        plot_box(w,box)
    
    Tk.mainloop()


def plot_box(canvas,box):
    print canvas.width,canvas.height,canvas.left,canvas.right,canvas.bottom,canvas.top
    assert(box.dimension()>=2)
    (xl,yl,xu,yu) = (box[0].lower(),box[1].lower(),box[0].upper(),box[1].upper())
    c=canvas
    xl=int((xl-c.left)/(c.right-c.left)*c.width)
    xu=int((xu-c.left)/(c.right-c.left)*c.width)
    yl=int((yl-c.top)/(c.bottom-c.top)*c.height)
    yu=int((yu-c.top)/(c.bottom-c.top)*c.height)
    print "  ",xl,yl,xu,yu
    canvas.create_rectangle(xl, yl, xu, yu, fill="blue")

   
if __name__=='__main__':
    boxes = [ Box([[0.1,0.2],[0.3,0.4]]), Box([[0.2,0.4],[0.4,0.7]]) ]
    plot(boxes)
