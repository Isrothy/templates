import java.math.BigInteger;
public class BigIntegerExample {
    public static void main(String[] args) {
        BigInteger one = BigInteger.ONE;
        BigInteger two = BigInteger.TWO;
        BigInteger zero = BigInteger.ZERO;
        BigInteger ten = BigInteger.TEN;
        BigInteger i114514 = BigInteger.valueOf(114514);
        BigInteger fromString = new BigInteger("19260817");
        BigInteger added = one.add(two);
        BigInteger subtracted = added.subtract(ten);
        BigInteger multiplied = subtracted.multiply(i114514);
        BigInteger divided = multiplied.divide(subtracted);
        BigInteger reminder = fromString.remainder(divided);
        System.out.println(reminder);
    }
}
