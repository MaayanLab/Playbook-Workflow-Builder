export class ResponseCodedError extends Error {
  constructor(public error_code: number, message?: string) {
    super(message)
    Object.setPrototypeOf(this, ResponseCodedError.prototype);
  }
}

/**
 * The method is not supported
 */
export class UnsupportedMethodError extends ResponseCodedError {
  constructor() {
    super(405, 'Unsupported method')
  }
}

/**
 * The resource is not found
 */
export class NotFoundError extends ResponseCodedError {
  constructor() {
    super(404, 'Not Found')
  }
}

/**
 * The resource is not found
 */
export class UnauthorizedError extends ResponseCodedError {
  constructor() {
    super(401, 'Unauthorized')
  }
}

/**
 * This timeout error is used to ensure we don't wait too long for dependencies
 *  fortunately, even if it occurs the job will requeue still making progress.
 */
export class TimeoutError extends ResponseCodedError {
  constructor() {
    super(504, 'Timeout reached')
  }
}

/**
 * This error occurs when the input node is not populated yet
 */
export class UnboundError extends ResponseCodedError {
  constructor() {
    super(422, 'Refusing to submit unbound variable')
  }
}
