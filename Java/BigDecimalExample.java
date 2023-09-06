import java.math.BigDecimal;
import java.math.BigInteger;
class BigIntegerExample {
    public static void main(String[] args) {
        BigDecimal one = BigDecimal.ONE;
        BigDecimal ten = BigDecimal.TEN;
        BigDecimal zero = BigDecimal.ZERO;
        BigDecimal fromInteger = new BigDecimal(BigInteger.valueOf(114514));
        BigDecimal fromString = new BigDecimal("1926.0817");
    }
}
