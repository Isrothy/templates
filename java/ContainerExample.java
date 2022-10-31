import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

class ContainerExample {
    public static void main(String[] args) {
        HashMap map = new HashMap<Integer, String>();
        map.put(10, "Ten");
        map.put(1, "One");
        System.out.println(map.containsKey(10));
        System.out.println(map.get(1));
        map.remove(10);

        HashSet set = new HashSet<Integer>();
        set.add(1);
        set.add(233);
        set.add(114514);
        set.remove(233);

        ArrayList list = new ArrayList<>();
        list.add(123);
        list.add(1, 10);
        System.out.println(list.get(1));
    }
}
